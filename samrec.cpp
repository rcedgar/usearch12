#include "myutils.h"
#include "samrec.h"
#include "seqinfo.h"
#include <limits>

/***
Tabbed fields
	1	QNAME	String [!-?A-~]f1,255g			Read label. Assumed sorted.
	2	FLAG	Int [0,216-1]					Bitwise FLAG.
	3	RNAME	String \*|[!-()+-<>-~][!-~]*	Target label.
	4	POS		Int [0,229-1]					1-based position in target.
	5	MAPQ	Int [0,28-1]					MAPping Quality. Q = -10 log10 Pe (error=wrong locus). Totally useless.
	6	CIGAR	String \*|([0-9]+[MIDNSHPX=])+	CIGAR string. Bowtie does not use = or X, only uses M.
	7	RNEXT	String \*|=|[!-()+-<>-~][!-~]*	Ref. name of the mate/next segment.
	8	PNEXT	Int [0,229-1]					Position of the mate/next segment.
	9	TLEN	Int [-229+1,229-1]				Template LENgth. Bowtie sets to 0 for unpaired reads.
	10	SEQ		String \*|[A-Za-z=.]+			Read segment SEQuence.
	11	QUAL	String [!-~]+					ASCII of Phred-scaled base QUALity+33.
	12,13...									Tags (tab-separated).

Bitwise flags:
	0x1		1		template having multiple segments in sequencing
	0x2		2		each segment properly aligned according to the aligner
	0x4		4		segment unmapped
	0x8		8		next segment in the template unmapped
	0x10	16		SEQ being reverse complemented
	0x20	32		SEQ of the next segment in the template being reversed
	0x40	64		the first segment in the template
	0x80	128		the last segment in the template
	0x100	256		secondary alignment (BWA uses 0x800 instead)
	0x200	512		not passing quality control
	0x400	1024	PCR or optical duplicate
	0x800	2048	secondary alignment (BWA only, not in SAM spec).

CIGAR
	MIDNSHPX=
	M 0 alignment match (can be a sequence match or mismatch)
	I 1 insertion to the reference
	D 2 deletion from the reference
	N 3 skipped region from the reference
	S 4 soft clipping (clipped sequences present in SEQ)
	H 5 hard clipping (clipped sequences NOT present in SEQ)
	P 6 padding (silent deletion from padded reference)
	= 7 sequence match
	X 8 sequence mismatch

Tags:

AS:i:-44	XS:i:-44	XN:i:0	XM:i:6	XO:i:1	XG:i:1	NM:i:7	MD:Z:25T5T34G51T0T0C0	YT:Z:UU

AS	Alignment score. Can be negative. Can be greater than 0 in --local mode (but not in
	--end-to-end mode).

XS	Alignment score for second-best alignment. Can be negative. Can be greater than 0 in 
	--local mode (but not in --end-to-end mode).

YS	Alignment score for opposite mate in the paired-end alignment.

XN	The number of ambiguous bases in the reference covering this alignment.

XM	The number of mismatches in the alignment.

NM	Edit distance to the reference, including ambiguous bases but excluding clipping.

XO	The number of gap opens, for both read and reference gaps, in the alignment.

MD	String for mismatching positions. Regex : [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*.

XG	The number of gap extensions, for both read and reference gaps, in the alignment. 

YT	UU=not part of a pair. CP=pair aligned concordantly. DP=part of a pair and the pair
	aligned discordantly. UP=pair failed to aligne either concordantly or discordantly.

MD	A string representation of the mismatched reference bases in the alignment.

The MD field aims to achieve SNP/indel calling without looking at the reference. For example,
a string `10A5^AC6' means from the leftmost reference base in the alignment, there are 10 matches
followed by an A on the reference which is different from the aligned read base; the next 5
reference bases are matches followed by a 2bp deletion from the reference; the deleted sequence
is AC; the last 6 bases are matches. The MD field ought to match the CIGAR string.

Reverse-complemented hits.
	For minus strand hit, POS is position of lowest target coordinate
	  in plus strand target sequence. This is the position aligned to the
	  first base in the rev-comp'd read / last base in the original read.
	The read sequence in the SAM file is revcomp'd.
	Conceptually, SAM rev-comps the read but not the target sequence.
***/

void SeqToFasta(FILE *f, const byte *Seq, unsigned L, const char *Label);

static double GetBlastnEvalue(unsigned QueryLength, double BlastnScore, double DBLength)
	{
// From gapped BLASTN:
	const double GAPPED_LAMBDA = 1.28;
	const double GAPPED_K = 0.46;

	static const double Log2 = log(2.0);

	double NM = double(QueryLength*DBLength);
	double BitScore = (BlastnScore*GAPPED_LAMBDA - GAPPED_K)/Log2;
	double Evalue = NM/pow(2, BitScore);

// Underflow tests
	if (Evalue != Evalue) // test for NaN
		return 0.0;

// Probably conservative, but better safe than sorry...
	if (numeric_limits<double>::has_infinity && Evalue == numeric_limits<double>::infinity())
		return 0.0;

	return Evalue;
	}

void CigarToVecs(const string &Cigar, vector<unsigned> &Ints,
  vector<char> &Chars)
	{
	Ints.clear();
	Chars.clear();

	const unsigned CL = SIZE(Cigar);
	if (CL == 0)
		return;

	unsigned Int = 0;
	asserta(isdigit(Cigar[0]) && !isdigit(Cigar[CL-1]));
	for (unsigned i = 0; i < CL; ++i)
		{
		char c = Cigar[i];
		if (isdigit(c))
			Int = Int*10 + (c - '0');
		else
			{
			asserta(i > 0 && isdigit(Cigar[i-1]));
			asserta(Int > 0);
			Ints.push_back(Int);
			Chars.push_back(c);
			Int = 0;
			}
		}
	}

// bowtie bug: bowtie2 produces bizarre CIGARs like this:
//	271S576M291I1M2686S
//          ^^^^^^
// Note: 291I 1M, very long gap followed by one or two matches.

SAMRec *SAMRec::MakeCopy() const
	{
	SAMRec *pr = new SAMRec;
	if (pr == 0)
		Die("SAMRec::MakeCopy, out of memory");

	pr->FromLine(m_Line.c_str());
	return pr;
	}

void SAMRec::ReadSeqToFasta(FILE *f) const
	{
	if (f == 0)
		return;
	asserta(m_ReadSeq != 0);

	SeqToFasta(f, (const byte *) m_ReadSeq, GetReadSeqLength(), m_ReadLabel);
	}

void SAMRec::ToBlast6(FILE *f)
	{
	if (f == 0)
		return;

	if (IsUnmapped())
		{
		if (opt(output_no_hits))
			{
			fprintf(f,
			  "%s"						//  1  Query label
			  "\t*"						//  2  Target label
			  "\t0"						//  3  %id
			  "\t0"						//  4  Alignment length (cols)
			  "\t0"						//  5  Mismatches
			  "\t0"						//  6  Opens
			  "\t0"						//  7  Qlo
			  "\t0"						//  8  Qhi
			  "\t0"						//  9  Tlo
			  "\t0"						// 10  Thi
			  "\t0"						// 10  Thi
			  "\t*"						// 11  Evalue
			  "\t0\n",					// 12  BitScore
			  m_ReadLabel);
			}
		return;
		}

	unsigned Pos = GetReadLo();
	unsigned ReadLo = Pos + 1;
	unsigned ReadHi = GetReadHi() + 1;
	unsigned TargetLo = GetTargetLo() + 1;
	unsigned TargetHi = GetTargetHi() + 1;
	double PctId = GetPctId();

	fprintf(f, "%s", m_ReadLabel);			// 1
	fprintf(f, "\t%s", m_TargetLabel);		// 2

	if (PctId == DBL_MAX)
		fprintf(f, "\t*");					// 3
	else
		fprintf(f, "\t%.1f", PctId);		// 3

	fprintf(f, "\t%u", m_ColCount);			// 4
	fprintf(f, "\t%u", m_GapOpens);			// 5
	fprintf(f, "\t%u", m_GapExts);			// 6

	fprintf(f, "\t%u", ReadLo);				// 7
	fprintf(f, "\t%u", ReadHi);				// 8

	if (m_RStrand)
		{
		fprintf(f, "\t%u", TargetLo);		// 9
		fprintf(f, "\t%u", TargetHi);		// 10
		}
	else
		{
		fprintf(f, "\t%u", TargetHi);		// 9
		fprintf(f, "\t%u", TargetLo);		// 10
		}

	if (m_Score == UINT_MAX)
		{
		fprintf(f, "\t*");					// 11
		fprintf(f, "\t*");					// 12
		}
	else
		{
		fprintf(f, "\t%.2g", GetEvalue());	// 11
		fprintf(f, "\t%u", m_Score);		// 12
		}

	fprintf(f, "\n");
	}

void SAMRec::ToFile(FILE *f) const
	{
	if (f == 0)
		return;

	unsigned FieldCount = SIZE(m_Fields);
	for (unsigned i = 0; i < FieldCount; ++i)
		{
		if (i > 0)
			fputc('\t', f);
		fputs(m_Fields[i].c_str(), f);
		}
	fputc('\n', f);
	}

unsigned SAMRec::GetReadLength() const
	{
	return CigarToReadLength(m_Cigar);
	}

double SAMRec::GetReadCovPct() const
	{
	return GetPct(m_Identities + m_Mismatches, GetReadLength());
	}

double SAMRec::GetEvalue() const
	{
	double BlastnScore = m_Identities - 2.0*(m_Mismatches + m_GapOpens + m_GapExts);
	unsigned L = GetReadLength();
	double Evalue = GetBlastnEvalue(L, BlastnScore, opt(ka_dbsize));
	return Evalue;
	}

unsigned SAMRec::CigarToColCount(const string &sCigar)
	{
	const char *Cigar = sCigar.c_str();
	unsigned Cols = 0;
	unsigned Int = 0;
	for (const char *p = Cigar; ; ++p)
		{
		char c = *p;
		if (c == 0)
			break;
		if (isdigit(c))
			Int = Int*10 + (c - '0');
		else if (isalpha(c) || c == '=')
			{
			asserta(Int > 0);
			if (c != 'S' && c != 'H')
				Cols += Int;
			Int = 0;
			}
		else
			Die("Unrecognized letter '%c' in CIGAR '%s'", c, Cigar);
		}
	return Cols;
	}

unsigned SAMRec::CigarToReadLo(const string &sCigar)
	{
	const char *Cigar = sCigar.c_str();
	unsigned Cols = 0;
	unsigned Int = 0;
	for (const char *p = Cigar; ; ++p)
		{
		char c = *p;
		if (c == 0)
			break;
		if (isdigit(c))
			Int = Int*10 + (c - '0');
		else if (isalpha(c) || c == '=')
			{
			if (c == 'S' || c == 'H')
				return Int;
			else
				return 0;
			Int = 0;
			}
		else
			Die("Unrecognized letter '%c' in CIGAR '%s'", c, Cigar);
		}
	return 0;
	}

unsigned SAMRec::CigarToReadSeqLength(const string &sCigar)
	{
	const char *Cigar = sCigar.c_str();
	unsigned ReadSeqLength = 0;
	unsigned Int = 0;
	for (const char *p = Cigar; ; ++p)
		{
		char c = *p;
		if (c == 0)
			break;
		if (isdigit(c))
			Int = Int*10 + (c - '0');
		else if (isalpha(c) || c == '=')
			{
			switch (c)
				{
			case 'D':
			case 'H':
				break;
			case 'S':
			case 'M':
			case 'I':
				ReadSeqLength += Int;
				break;

			default:
				Die("Unrecognized letter '%c' in CIGAR '%s'", c, Cigar);
				}
			Int = 0;
			}
		}
	return ReadSeqLength;
	}

unsigned SAMRec::CigarToReadLength(const string &sCigar)
	{
	const char *Cigar = sCigar.c_str();
	unsigned ReadLength = 0;
	unsigned Int = 0;
	for (const char *p = Cigar; ; ++p)
		{
		char c = *p;
		if (c == 0)
			break;
		if (isdigit(c))
			Int = Int*10 + (c - '0');
		else if (isalpha(c) || c == '=')
			{
			switch (c)
				{
			case 'H':
			case 'S':
			case 'M':
			case 'I':
				ReadLength += Int;
				break;

			case 'D':
				break;

			default:
				Die("Unrecognized letter '%c' in CIGAR '%s'", c, Cigar);
				}
			Int = 0;
			}
		}
	return ReadLength;
	}

unsigned SAMRec::CigarToReadSeqLo(const string &sCigar)
	{
	const char *Cigar = sCigar.c_str();
	unsigned Cols = 0;
	unsigned Int = 0;
	for (const char *p = Cigar; ; ++p)
		{
		char c = *p;
		if (c == 0)
			break;
		if (isdigit(c))
			Int = Int*10 + (c - '0');
		else if (isalpha(c) || c == '=')
			{
			if (c == 'S')
				return Int;
			else
				return 0;
			Int = 0;
			}
		else
			Die("Unrecognized letter '%c' in CIGAR '%s'", c, Cigar);
		}
	return 0;
	}

void SAMRec::FromLine(const char *Line)
	{
	FromLine(string(Line));
	}

void SAMRec::Clear()
	{
	m_Line.clear();
	m_Cigar.clear();
	m_CigarChars.clear();
	m_CigarInts.clear();

	m_RRow.Size = 0;
	m_ARow.Size = 0;
	m_TRow.Size = 0;
	m_TSeg.Size = 0;

	m_TargetLabel = 0;
	m_ReadSeq = 0;
	m_TargetLabel = 0;
	m_Score = UINT_MAX;
	m_BitFlags = 0;
	m_TargetLo = UINT_MAX;
	m_ReadSeqLo = UINT_MAX;
	m_ReadLo = UINT_MAX;
	m_ReadSeqLength = UINT_MAX;
	m_Identities = UINT_MAX;
	m_Mismatches = UINT_MAX;
	m_Gaps = UINT_MAX;
	m_GapOpens = UINT_MAX;
	m_GapExts = UINT_MAX;
	m_ColCount = UINT_MAX;
	m_Qual = 0;
	m_MD = 0;
	m_Evalue = DBL_MAX;
	m_RStrand = false;
	m_IsHardClipped = false;
	}

void SAMRec::FromLine(const string &Line)
	{
	Clear();

	m_Line = Line;

	Split(m_Line, m_Fields, '\t');
	unsigned FieldCount = SIZE(m_Fields);
	if (FieldCount < 11)
		Die("Invalid SAM record, got %u tabbed fields must be >= 11", FieldCount);

	m_ReadLabel = m_Fields[0].c_str();
	m_BitFlags = atou(m_Fields[1]);

// yara bug -- unmap if CIGAR is *
	if (m_Fields[5] == "*")
		m_BitFlags = 4;

	m_ReadSeq = m_Fields[9].c_str();
	m_ReadSeqLength = ustrlen(m_ReadSeq);
	m_Mapq = UINT_MAX;

	if (IsUnmapped())
		return;

	m_Mapq = atou(m_Fields[4]);
	if (m_Mapq == 255)
		m_Mapq = UINT_MAX;

	m_TargetLabel = m_Fields[2].c_str();

	m_TargetLo = atou(m_Fields[3]);
	asserta(m_TargetLo > 0);
	--m_TargetLo; // convert to 0-based
	m_Cigar = m_Fields[5];
	m_Qual = m_Fields[10].c_str();
	m_MD = 0;
	m_Identities = UINT_MAX;
	m_Mismatches = UINT_MAX;
	m_GapOpens = UINT_MAX;
	m_Gaps = UINT_MAX;
	m_Score = UINT_MAX;
	m_RStrand = ((m_BitFlags & 16) == 0);

	for (unsigned i = 11; i < FieldCount; ++i)
		{
		const char *f = m_Fields[i].c_str();
		if (StartsWith(f, "MD:Z:"))
			m_MD = f + 5;
		else if (StartsWith(f, "XM:i:"))
			m_Mismatches = atou(f + 5);
		else if (StartsWith(f, "XO:i:"))
			m_GapOpens = atou(f + 5);
		else if (StartsWith(f, "XG:i:"))
			m_Gaps = atou(f + 5);
		else if (StartsWith(f, "AS:i:"))
			m_Score = atou(f + 5);
		}

	m_ColCount = CigarToColCount(m_Cigar);

	if (m_GapOpens == UINT_MAX)
		m_GapOpens = CalcGapOpens();
	if (m_Gaps == UINT_MAX)
		m_Gaps = CalcGaps();

	if (m_GapOpens != UINT_MAX && m_Gaps != UINT_MAX)
		{
		asserta(m_Gaps >= m_GapOpens);
		m_GapExts = m_Gaps - m_GapOpens;
		}

	unsigned NonIdentities = UINT_MAX;
	if (m_Mismatches != UINT_MAX && m_Gaps != UINT_MAX)
		{
		NonIdentities = m_Mismatches + m_Gaps;
		m_Identities = m_ColCount - NonIdentities;
		}

	m_ReadLo = CigarToReadLo(m_Cigar);
	m_ReadSeqLo = CigarToReadSeqLo(m_Cigar);
	m_IsHardClipped = (strchr(m_Cigar.c_str(), 'H') != 0);
	}

const char *SAMRec::GetRRow()
	{
	if (m_RRow.Size > 0)
		return m_RRow.Data;

	m_RRow.Alloc(m_ColCount + 16);
	char *RRow = m_RRow.Data;

	const unsigned ReadSeqLength = GetReadSeqLength();
	asserta(m_ReadSeq[ReadSeqLength] == 0);

	const unsigned CL = SIZE(m_Cigar);
	unsigned Int = 0;
	unsigned Col = 0;
	bool FirstS = true;
	unsigned ReadSeqPos = 0;
	for (unsigned i = 0; i < CL; ++i)
		{
		char c = m_Cigar[i];
		if (isdigit(c))
			Int = Int*10 + (c - '0');
		else
			{
			switch (c)
				{
			case 'H':
				break;

			case 'S':
				if (FirstS)
					{
					ReadSeqPos = Int;
					FirstS = false;
					}
				else
					ReadSeqPos += Int;
				break;

			case 'M':
			case 'I':
				for (unsigned j = 0; j < Int; ++j)
					{
					if (ReadSeqPos >= ReadSeqLength)
						{
						Log("\n");
						Log("%s\n", m_Line.c_str());
						Die("ReadPos %u >= ReadSeqLength %u", ReadSeqPos, ReadSeqLength);
						}
					RRow[Col++] = toupper(m_ReadSeq[ReadSeqPos++]);
					}
				break;

			case 'D':
				for (unsigned j = 0; j < Int; ++j)
					RRow[Col++] = '-';
				break;

			default:
				Die("Unexpected char '%c' in CIGAR '%s'", c, m_Cigar.c_str());
				}
			Int = 0;
			}
		}
	asserta(Col == m_ColCount);
	RRow[Col] = 0;
	m_RRow.Size = m_ColCount;
	return RRow;
	}

const char *SAMRec::GetTRow()
	{
	if (m_MD == 0)
		Die("SAMRec::GetTSeg(), SAM record missing MD, cannot convert to alignment");

	if (m_TSeg.Size > 0)
		{
		const char *TSeg = m_TSeg.Data;
		asserta(TSeg != 0);
		return TSeg;
		}

	m_TSeg.Alloc(m_ColCount+16);
	char *TSeg = m_TSeg.Data;
	memset(TSeg, '?', m_ColCount);

	unsigned Int = 0;
	const unsigned L = ustrlen(m_MD);
	bool Deleting = false;
	unsigned Pos = 0;
	for (unsigned i = 0; i < L; ++i)
		{
		char c = m_MD[i];
		if (isdigit(c))
			{
			Deleting = false;
			Int = Int*10 + (c - '0');
			}
		else if (isalpha(c))
			{
			for (unsigned j = 0; j < Int; ++j)
				TSeg[Pos++] = '.';

			Int = 0;
			TSeg[Pos++] = c;
			}
		else
			{
			asserta(c == '^');
			Deleting = true;
			}
		}

	for (unsigned j = 0; j < Int; ++j)
		TSeg[Pos++] = '.';

	TSeg[Pos] = 0;

	const unsigned ReadSeqLength = GetReadSeqLength();
	const unsigned CL = SIZE(m_Cigar);
	Int = 0;
	bool FirstS = true;
	unsigned ReadSeqPos = 0;
	unsigned TargetSeqPos = 0;
	for (unsigned i = 0; i < CL; ++i)
		{
		char c = m_Cigar[i];
		if (isdigit(c))
			Int = Int*10 + (c - '0');
		else
			{
			switch (c)
				{
			case 'H':
				break;

			case 'S':
				if (FirstS)
					{
					ReadSeqPos = Int;
					FirstS = false;
					}
				else
					ReadSeqPos += Int;
				break;

			case 'M':
				for (unsigned j = 0; j < Int; ++j)
					{
					if (ReadSeqPos >= ReadSeqLength)
						{
						Log("\n");
						Log("%s\n", m_Line.c_str());
						Die("ReadPos %u >= ReadSeqLength %u", ReadSeqPos, ReadSeqLength);
						}

					if (TSeg[TargetSeqPos] == '.')
						TSeg[TargetSeqPos] = m_ReadSeq[ReadSeqPos];

					++TargetSeqPos;
					++ReadSeqPos;
					}
				break;

			case 'I':
				ReadSeqPos += Int;
				break;

			case 'D':
				for (unsigned j = 0; j < Int; ++j)
					{
					char c = TSeg[TargetSeqPos++];
					asserta(c != '.' && c != '?');
					}
				break;

			default:
				Die("Unexpected char '%c' in CIGAR '%s'", c, m_Cigar.c_str());
				}
			Int = 0;
			}
		}

	m_TRow.Alloc(m_ColCount+16);
	char *TRow = m_TRow.Data;
	unsigned Col = 0;
	TargetSeqPos = 0;
	for (unsigned i = 0; i < CL; ++i)
		{
		char c = m_Cigar[i];
		if (isdigit(c))
			Int = Int*10 + (c - '0');
		else
			{
			switch (c)
				{
			case 'H':
			case 'S':
				break;

			case 'M':
			case 'D':
				for (unsigned j = 0; j < Int; ++j)
					TRow[Col++] = TSeg[TargetSeqPos++];
				break;

			case 'I':
				for (unsigned j = 0; j < Int; ++j)
					TRow[Col++] = '-';
				break;

			default:
				Die("Unexpected char '%c' in CIGAR '%s'", c, m_Cigar.c_str());
				}
			Int = 0;
			}
		}
	TRow[Col] = 0;
	if (Col != m_ColCount)
		{
		Log("RRow %s\n", GetRRow());
		Log("TRow %s\n", TRow);
		Die("Col %u, m_ColCount %u\n", Col, m_ColCount);
		}
	return TRow;
	}

const char *SAMRec::GetARow()
	{
	if (m_ARow.Size > 0)
		{
		const char *ARow = m_ARow.Data;
		asserta(ARow != 0);
		asserta(ARow[m_ColCount] == 0);
		return ARow;
		}

	asserta(m_ColCount > 0 && m_ColCount != UINT_MAX);
	const char *RRow = GetRRow();
	const char *TRow = GetTRow();
	asserta(RRow[m_ColCount] == 0);
	asserta(TRow[m_ColCount] == 0);

	m_ARow.Alloc(m_ColCount+16);
	char *ARow = m_ARow.Data;

	for (unsigned Col = 0; Col < m_ColCount; ++Col)
		{
		if (RRow[Col] == TRow[Col])
			ARow[Col] = '|';
		else
			ARow[Col] = ' ';
		}
	ARow[m_ColCount] = 0;
	m_ARow.Size = m_ColCount;
	return ARow;
	}

unsigned SAMRec::GetInsertCount()
	{
	SetCigarVecs();
	const unsigned N = SIZE(m_CigarChars);
	unsigned n = 0;
	for (unsigned i = 0; i < N; ++i)
		if (m_CigarChars[i] == 'I')
			n += m_CigarInts[i];
	return n;
	}

unsigned SAMRec::GetDeleteCount()
	{
	SetCigarVecs();
	const unsigned N = SIZE(m_CigarChars);
	unsigned n = 0;
	for (unsigned i = 0; i < N; ++i)
		if (m_CigarChars[i] == 'D')
			n += m_CigarInts[i];
	return n;
	}

unsigned SAMRec::GetTargetHi()
	{
	return m_TargetLo + GetTargetSegLength() - 1;
	}

unsigned SAMRec::GetReadSegLength()
	{
	unsigned DeleteCount = GetDeleteCount();
	asserta(DeleteCount < m_ColCount);
	unsigned L = m_ColCount - DeleteCount;
	return L;
	}

unsigned SAMRec::GetTargetSegLength()
	{
	unsigned InsertCount = GetInsertCount();
	asserta(InsertCount < m_ColCount);
	unsigned L = m_ColCount - InsertCount;
	return L;
	}

void SAMRec::PrAln(FILE *f)
	{
	if (f == 0)
		return;

	const char *RRow = GetRRow();
	const char *TRow = GetTRow();
	const char *ARow = GetARow();

	unsigned Identities = 0;
	unsigned ReadAlignedBaseCount = 0;
	for (unsigned Col = 0; Col < m_ColCount; ++Col)
		{
		char r = RRow[Col];
		char t = TRow[Col];
		if (r != '-')
			{
			++ReadAlignedBaseCount;
			if (r == t)
				++Identities;
			}
		}

	const unsigned COLS = 80;
	unsigned Rows = (m_ColCount + COLS - 1)/COLS;

	unsigned QStart = (m_RStrand ? GetReadLo() : GetReadHi()) + 1;
	unsigned QEnd = (m_RStrand ? GetReadHi() : GetReadLo()) + 1;
	unsigned QMax = max(QStart, QEnd);

	unsigned TStart = GetTargetLo() + 1;
	unsigned TEnd = GetTargetHi() + 1;

	unsigned MaxPos = max(QMax, TEnd);

	char Tmp[16];
	sprintf(Tmp, "%u", MaxPos);
	unsigned ndig = (unsigned) strlen(Tmp);

	Pr(f, "\n");

	Pr(f, "Query  >%s, length %u nt\n", m_ReadLabel, GetReadLength());
	Pr(f, "Target >%s", m_TargetLabel);
	Pr(f, "\n");

	unsigned QPos = QStart;
	unsigned TPos = TStart;
	for (unsigned Row = 0; Row < Rows; ++Row)
		{
		unsigned ColLo = Row*COLS;
		unsigned n = m_ColCount - ColLo;
		if (n > COLS)
			n = COLS;

		Pr(f, "\n");
		Pr(f, "Qry %*u %c ", ndig, QPos, pom(m_RStrand));
		WriteStdioFile(f, RRow + ColLo, n);
		for (unsigned i = 0; i < n; ++i)
			{
			unsigned Col = ColLo + i;
			char q = RRow[Col];
			if (q != '-')
				{
				if (m_RStrand)
					++QPos;
				else
					--QPos;
				}
			}
		Pr(f, " %u\n", m_RStrand ? QPos - 1 : QPos + 1);

		Pr(f, "      %*.*s ", ndig, ndig, "");
		WriteStdioFile(f, ARow + ColLo, n);
		Pr(f, "\n");

		Pr(f, "Tgt %*u + ", ndig, TPos);
		WriteStdioFile(f, TRow + ColLo, n);
		for (unsigned i = 0; i < n; ++i)
			{
			unsigned Col = ColLo + i;
			char t = TRow[Col];
			if (t != '-')
				++TPos;
			}
		Pr(f, " %u\n", TPos - 1);
		}

	unsigned Cols = GetColCount();
	unsigned Gaps = m_GapOpens + m_GapExts;

	double ReadCovPct = GetPct(ReadAlignedBaseCount, m_ReadSeqLength);
	double PctId = GetPct(Identities, m_ColCount);

	Pr(f, "\n");
	Pr(f, "MAPQ %u", m_Mapq);
	Pr(f, ", score %u", m_Score);
	Pr(f, ", cols %u", Cols);
	Pr(f, ", cov %.1f%%", ReadCovPct);
	Pr(f, ", ids %u (%.1f%%)", Identities, PctId);
	Pr(f, ", gaps %u (%.1f%%)", Gaps, GetPct(Gaps, Cols));
	if (optset_ka_dbsize)
		Pr(f, ", E-value %.2g", GetEvalue());
	Pr(f, "\n");
	}

double SAMRec::GetPctId2()
	{
	const char *RRow = GetRRow();
	const char *TRow = GetTRow();
	unsigned Cols = 0;
	unsigned n = 0;
	for (unsigned i = 0; i < m_ColCount; ++i)
		{
		if (RRow[i] == 'N' || RRow[i] == '-' || TRow[i] == '-')
			continue;
		++Cols;
		if (RRow[i] == TRow[i])
			++n;
		}
	return GetPct(n, Cols);
	}

void SAMRec::GetColToTargetPos(vector<unsigned> &ColToTargetPos)
	{
	ColToTargetPos.clear();

	const char *RRow = GetRRow();
	const char *TRow = GetTRow();

	unsigned TargetPos = m_TargetLo;
	unsigned MaxTargetPos = TargetPos;
	for (unsigned Col = 0; Col < m_ColCount; ++Col)
		{
		char t = TRow[Col];

		if (t != '-')
			{
		 // convert to 0-based
			ColToTargetPos.push_back(TargetPos);
			MaxTargetPos = TargetPos;
			++TargetPos;
			}
		else
			ColToTargetPos.push_back(UINT_MAX);
		}
	unsigned TargetHi = GetTargetHi();
	asserta(MaxTargetPos == TargetHi);
	asserta(SIZE(ColToTargetPos) == m_ColCount);
	}

void SAMRec::GetPath(string &Path, bool ReadIsQuery)
	{
	Path.clear();
	const char *RRow = GetRRow();
	const char *TRow = GetTRow();
	unsigned ColCount = GetColCount();

	char ReadGapSymbol = (ReadIsQuery ? 'D' : 'I');
	char TargetGapSymbol = (ReadIsQuery ? 'I' : 'D');

	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		char r = RRow[Col];
		char t = TRow[Col];

		if (isalpha(t) && isalpha(t))
			Path.push_back('M');
		else if (r == '-')
			Path.push_back(ReadGapSymbol);
		else if (t == '-')
			Path.push_back(TargetGapSymbol);
		else
			asserta(false);
		}
	}

void SAMRec::ValidateReadSeqLength()
	{
	unsigned L1 = ustrlen(m_ReadSeq);
	unsigned L2 = CigarToReadSeqLength(m_Cigar);
	asserta(L1 == L2);
	}

void SAMRec::ValidateTRow(SeqInfo *Target)
	{
	const byte *T = Target->m_Seq;
	unsigned TL = Target->m_L;
	unsigned TLo = GetTargetLo();
	unsigned THi = GetTargetHi();
	if (THi < TLo || THi >= TL)
		Die("Bad rec THi %u, TLo %u, TL %u >%s", THi, TLo, TL, Target->m_Label);
	asserta(THi < TL);
	const char *TRow = GetTRow();
	unsigned TPos = TLo;
	unsigned TSegLength = GetTargetSegLength();
	unsigned ColCount = m_ColCount;
	for (unsigned Col = 0; Col < ColCount; ++Col)
		{
		char c = TRow[Col];
		if (c == '-')
			continue;
		char t = T[TPos++];
		if (toupper(t) != toupper(c))
			{
			Log("\n");
			Log("%s\n", m_Line.c_str());
			Log("RRow %*.*s\n", ColCount, ColCount, GetRRow());
			Log("TRow %*.*s\n", ColCount, ColCount, TRow);
			Log("TSeg %*.*s\n", TSegLength, TSegLength, T + TLo);
			Warning("Mismatch: Col %u, c %c, t %c\n", Col, c, t);
			}
		}
	}

unsigned SAMRec::CalcGapOpens()
	{
	SetCigarVecs();
	const unsigned N = SIZE(m_CigarChars);
	unsigned n = 0;
	for (unsigned i = 0; i < N; ++i)
		{
		char c = m_CigarChars[i];
		if (c == 'D' || c == 'I')
			++n;
		}
	return n;
	}

unsigned SAMRec::CalcGaps()
	{
	SetCigarVecs();
	const unsigned N = SIZE(m_CigarChars);
	unsigned n = 0;
	for (unsigned i = 0; i < N; ++i)
		{
		char c = m_CigarChars[i];
		if (c == 'D' || c == 'I')
			n += m_CigarInts[i];
		}
	return n;
	}
