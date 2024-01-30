#include "myutils.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "seqsource.h"
#include "seqdb.h"
#include "label.h"

static unsigned GetBarcodeIndex(const SeqDB &Barcodes, const byte *Q, unsigned QL,
  bool PrefixMatch, unsigned &MinDiffs)
	{
	const unsigned BarcodeCount = Barcodes.GetSeqCount();
	unsigned BestIndex = UINT_MAX;
	MinDiffs = UINT_MAX;
	for (unsigned BarcodeIndex = 0; BarcodeIndex < BarcodeCount; ++BarcodeIndex)
		{
		const byte *T = Barcodes.GetSeq(BarcodeIndex);
		unsigned TL = Barcodes.GetSeqLength(BarcodeIndex);
		if (QL < TL)
			continue;
		if (!PrefixMatch && TL != QL)
			continue;
		unsigned Diffs = 0;
		for (unsigned i = 0; i < TL; ++i)
			{
			byte q = toupper(Q[i]);
			byte t = toupper(T[i]);
			if (q != t)
				++Diffs;
			}
		if (Diffs < MinDiffs)
			{
			BestIndex = BarcodeIndex;
			MinDiffs = Diffs;
			}
		}
	return BestIndex;
	}

static void DemuxStripDual()
	{
	asserta(!optset_reverse);
	asserta(!optset_index);

	FILE *fTab = 0;
	FILE *fFa = 0;
	FILE *fFq = 0;

	if (optset_tabbedout)
		fTab = CreateStdioFile(opt(tabbedout));
	if (optset_fastaout)
		fFa = CreateStdioFile(opt(fastaout));
	if (optset_fastqout)
		fFq = CreateStdioFile(opt(fastqout));

	SeqSource *SS = MakeSeqSource(opt(fastx_demux));

	unsigned MaxDiffs = 0;
	if (optset_maxdiffs)
		MaxDiffs = opt(maxdiffs);

	SeqDB Barcodes;
	Barcodes.FromFasta(opt(barcodes));

	string Delim;
	if (optset_sample_delim)
		Delim = opt(sample_delim);

	const unsigned BarcodeCount = Barcodes.GetSeqCount();

	SeqInfo *SI = ObjMgr::GetSeqInfo();
	SeqInfo *SIRC = ObjMgr::GetSeqInfo();

	unsigned InCount = 0;
	unsigned OutCount = 0;
	ProgressStep(0, 1000, "Processing");
	bool HasQual = false;
	for (;;)
		{
		ProgressStep(SS->GetPctDoneX10(), 1000, "Demuxed %u / %u (%.1f%%)",
		  OutCount, InCount, GetPct(OutCount, InCount));

		bool Ok = SS->GetNext(SI);
		if (!Ok)
			break;

		string ReadLabel = SI->m_Label;
		++InCount;

		SI->GetRevComp(SIRC);

		unsigned Diffs;
		unsigned DiffsRC;
		unsigned BarcodeIndex = GetBarcodeIndex(Barcodes, SI->m_Seq, SI->m_L, true, Diffs);
		unsigned BarcodeIndexRC = GetBarcodeIndex(Barcodes, SIRC->m_Seq, SIRC->m_L, true, DiffsRC);

		string BarcodeLabel = "*";
		string BarcodeLabelRC = "*";
		if (Diffs <= MaxDiffs && DiffsRC <= MaxDiffs)
			{
			++OutCount;

			asserta(BarcodeIndex < BarcodeCount);
			asserta(BarcodeIndexRC < BarcodeCount);
			BarcodeLabel = Barcodes.GetLabel(BarcodeIndex);
			BarcodeLabelRC = Barcodes.GetLabel(BarcodeIndexRC);
			unsigned BarcodeLength = Barcodes.GetSeqLength(BarcodeIndex);
			unsigned BarcodeLengthRC = Barcodes.GetSeqLength(BarcodeIndexRC);
			asserta(BarcodeLength + BarcodeLengthRC < SI->m_L);

			string Acc;
			GetAccFromLabel(SI->m_Label, Acc);

			string NewLabel;
			string BC = BarcodeLabel + "." + BarcodeLabelRC;
			if (optset_sample_delim)
				{
				char Tmp[16];
				sprintf(Tmp, "%u", OutCount);
				NewLabel = BC + Delim + string(Tmp);
				}
			else
				NewLabel = Acc + string(";sample=") + BC + ";";

			SI->m_Label = NewLabel.c_str();
			SI->StripLeft(BarcodeLength);
			SI->StripRight(BarcodeLengthRC);
			SI->ToFastq(fFq);
			SI->ToFasta(fFa);
			}

		if (fTab != 0)
			{
			fprintf(fTab, "%s", SI->m_Label);
			fprintf(fTab, "\t%s", BarcodeLabel.c_str());
			fprintf(fTab, "\t%s", BarcodeLabelRC.c_str());
			fprintf(fTab, "\t%u", Diffs);
			fprintf(fTab, "\t%u", DiffsRC);
			fprintf(fTab, "\n");
			}
		}

	ProgressStep(999, 1000, "Demuxed %u / %u (%.1f%%)",
	  OutCount, InCount, GetPct(OutCount, InCount));

	CloseStdioFile(fTab);
	CloseStdioFile(fFa);
	CloseStdioFile(fFq);
	}

static void DemuxStrip()
	{
	asserta(!optset_reverse);
	asserta(!optset_index);

	FILE *fTab = 0;
	FILE *fFa = 0;
	FILE *fFq = 0;

	if (optset_tabbedout)
		fTab = CreateStdioFile(opt(tabbedout));
	if (optset_fastaout)
		fFa = CreateStdioFile(opt(fastaout));
	if (optset_fastqout)
		fFq = CreateStdioFile(opt(fastqout));

	SeqSource *SS = MakeSeqSource(opt(fastx_demux));

	unsigned MaxDiffs = 0;
	if (optset_maxdiffs)
		MaxDiffs = opt(maxdiffs);

	SeqDB Barcodes;
	Barcodes.FromFasta(opt(barcodes));

	string Delim;
	if (optset_sample_delim)
		Delim = opt(sample_delim);

	const unsigned BarcodeCount = Barcodes.GetSeqCount();

	SeqInfo *SI = ObjMgr::GetSeqInfo();

	unsigned InCount = 0;
	unsigned OutCount = 0;
	ProgressStep(0, 1000, "Processing");
	bool HasQual = false;
	for (;;)
		{
		ProgressStep(SS->GetPctDoneX10(), 1000, "Demuxed %u / %u (%.1f%%)",
		  OutCount, InCount, GetPct(OutCount, InCount));

		bool Ok = SS->GetNext(SI);
		if (!Ok)
			break;

		string ReadLabel = SI->m_Label;

		++InCount;

		unsigned Diffs;
		unsigned BarcodeIndex = GetBarcodeIndex(Barcodes, SI->m_Seq, SI->m_L, true, Diffs);
		if (Diffs >= MaxDiffs)
			{
			SI->RevCompInPlace();
			BarcodeIndex = GetBarcodeIndex(Barcodes, SI->m_Seq, SI->m_L, true, Diffs);
			}

		string BarcodeLabel = "*";
		if (Diffs <= MaxDiffs)
			{
			++OutCount;

			asserta(BarcodeIndex < BarcodeCount);
			BarcodeLabel = Barcodes.GetLabel(BarcodeIndex);
			unsigned BarcodeLength = Barcodes.GetSeqLength(BarcodeIndex);
			asserta(BarcodeLength < SI->m_L);

			string Acc;
			GetAccFromLabel(SI->m_Label, Acc);

			string NewLabel;
			if (optset_sample_delim)
				{
				char Tmp[16];
				sprintf(Tmp, "%u", OutCount);
				NewLabel = BarcodeLabel + Delim + string(Tmp);
				}
			else
				NewLabel = Acc + string(";sample=") + BarcodeLabel + ";";

			SI->StripLeft(BarcodeLength);
			SI->ToFastq(fFq);
			SI->ToFasta(fFa);
			}

		if (fTab != 0)
			{
			fprintf(fTab, "%s", SI->m_Label);
			fprintf(fTab, "\t%s", BarcodeLabel.c_str());
			fprintf(fTab, "\t%u", Diffs);
			fprintf(fTab, "\n");
			}
		}

	ProgressStep(999, 1000, "Demuxed %u / %u (%.1f%%)",
	  OutCount, InCount, GetPct(OutCount, InCount));

	CloseStdioFile(fTab);
	CloseStdioFile(fFa);
	CloseStdioFile(fFq);
	}

void cmd_fastx_demux()
	{
	if (!optset_barcodes)
		Die("missing -barcodes");

	if (opt(strip))
		{
		DemuxStrip();
		return;
		}
	else if (opt(stripd))
		{
		DemuxStripDual();
		return;
		}

	FILE *fTab = 0;
	FILE *fFa = 0;
	FILE *fFq = 0;
	FILE *fFq2 = 0;

	if (optset_tabbedout)
		fTab = CreateStdioFile(opt(tabbedout));
	if (optset_fastaout)
		fFa = CreateStdioFile(opt(fastaout));
	if (optset_fastqout)
		fFq = CreateStdioFile(opt(fastqout));
	if (optset_output2)
		fFq2 = CreateStdioFile(opt(output2));

	SeqSource *SS = MakeSeqSource(opt(fastx_demux));
	SeqSource *SS2 = 0;
	if (optset_reverse)
		SS2 = MakeSeqSource(opt(reverse));
	SeqSource *SSI = 0;
	if (optset_index)
		SSI = MakeSeqSource(opt(index));

	unsigned MaxDiffs = 0;
	if (optset_maxdiffs)
		MaxDiffs = opt(maxdiffs);

	SeqDB Barcodes;
	Barcodes.FromFasta(opt(barcodes));

	string Delim;
	if (optset_sample_delim)
		Delim = opt(sample_delim);

	const unsigned BarcodeCount = Barcodes.GetSeqCount();

	SeqInfo *SI = ObjMgr::GetSeqInfo();
	SeqInfo *SI2 = ObjMgr::GetSeqInfo();
	SeqInfo *SII = ObjMgr::GetSeqInfo();

	unsigned InCount = 0;
	unsigned OutCount = 0;
	ProgressStep(0, 1000, "Processing");
	bool HasQual = false;
	for (;;)
		{
		ProgressStep(SS->GetPctDoneX10(), 1000, "Demuxed %u / %u (%.1f%%)",
		  OutCount, InCount, GetPct(OutCount, InCount));

		bool Ok = SS->GetNext(SI);
		if (!Ok)
			break;
		if (SS2 != 0)
			{
			bool Ok2 = SS2->GetNext(SI2);
			if (!Ok2)
				Die("Premature end-of-file in reverse reads");

			if (!IlluminaLabelPairMatch(SI->m_Label, SI2->m_Label))
				{
				Log("\n");
				Log("%s\n", SI->m_Label);
				Log("%s\n", SI2->m_Label);
				Die("Forward/reverse label mismatch");
				}
			}

		string ReadLabel = SI->m_Label;
		if (SSI != 0)
			{
			bool OkI = SSI->GetNext(SII);
			if (!OkI)
				Die("Premature end-of-file in index reads");

			vector<string> Fields;
			vector<string> FieldsI;
			string IndexLabel = SII->m_Label;
			Split(ReadLabel, Fields);
			Split(IndexLabel, FieldsI);
			if (Fields[0] != FieldsI[0])
				Die("Label mismatch read=%s, index=%s", Fields[0].c_str(), FieldsI[0].c_str());
			}

		++InCount;

		unsigned Diffs;
		unsigned BarcodeIndex = UINT_MAX;
		if (SSI != 0)
			BarcodeIndex = GetBarcodeIndex(Barcodes, SII->m_Seq, SII->m_L, false, Diffs);
		else
			BarcodeIndex = GetBarcodeIndex(Barcodes, SI->m_Seq, SI->m_L, true, Diffs);

		string BarcodeLabel = "*";
		if (Diffs <= MaxDiffs)
			{
			++OutCount;

			asserta(BarcodeIndex < BarcodeCount);
			BarcodeLabel = Barcodes.GetLabel(BarcodeIndex);

			string Acc;
			GetAccFromLabel(SI->m_Label, Acc);

			string NewLabel;
			if (optset_sample_delim)
				{
				char Tmp[16];
				sprintf(Tmp, "%u", OutCount);
				NewLabel = BarcodeLabel + Delim + string(Tmp);
				}
			else
				NewLabel = Acc + string(";sample=") + BarcodeLabel + ";";

			SeqToFasta(fFa, SI->m_Seq, SI->m_L, NewLabel.c_str());
			SeqToFastq(fFq, SI->m_Seq, SI->m_L, SI->m_Qual, NewLabel.c_str());
			if (SI2 != 0)
				SeqToFastq(fFq2, SI2->m_Seq, SI2->m_L, SI2->m_Qual, NewLabel.c_str());
			}

		if (fTab != 0)
			{
			fprintf(fTab, "%s", SI->m_Label);
			fprintf(fTab, "\t%s", BarcodeLabel.c_str());
			fprintf(fTab, "\t%u", Diffs);
			fprintf(fTab, "\n");
			}
		}

	ProgressStep(999, 1000, "Demuxed %u / %u (%.1f%%)",
	  OutCount, InCount, GetPct(OutCount, InCount));

	CloseStdioFile(fTab);
	CloseStdioFile(fFa);
	CloseStdioFile(fFq);
	CloseStdioFile(fFq2);
	}
