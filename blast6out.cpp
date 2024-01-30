#include "myutils.h"
#include "seqinfo.h"
#include "alignresult.h"
#include "outputsink.h"
#include "pathinfo.h"

/**
Blast -outfmt 6:

	qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

Fields:
            qseqid means Query Seq-id 1
            sseqid means Subject Seq-id 2
            pident means Percentage of identical matches 3
            length means Alignment length 4
          mismatch means Number of mismatches 5
           gapopen means Number of gap openings 6
            qstart means Start of alignment in query 7
              qend means End of alignment in query 8
            sstart means Start of alignment in subject 9
              send means End of alignment in subject 10
            evalue means Expect value 11
          bitscore means Bit score 12
**/

void OutputSink::OutputBlast6(AlignResult *AR)
	{
	if (m_fBlast6 == 0)
		return;

	const char *QLabel = AR->GetQueryLabel();
	const char *TLabel = AR->GetTargetLabel();

	unsigned QLo = AR->GetQLo6();
	unsigned QHi = AR->GetQHi6();

	unsigned TLo = AR->GetTLo6();
	unsigned THi = AR->GetTHi6();

	unsigned AlnLength = AR->GetAlnLength();
	unsigned MismatchCount = AR->GetMismatchCount();
	unsigned OpenCount = AR->GetGapOpenCount();
	double PctId = AR->GetPctId();

	fprintf(m_fBlast6,
	  "%s"						//  1  Query label
	  "\t%s"					//  2  Target label
	  "\t%.1f"					//  3  %id
	  "\t%u"					//  4  Alignment length (cols)
	  "\t%u"					//  5  Mismatches
	  "\t%u"					//  6  Opens
	  "\t%u"					//  7  Qlo
	  "\t%u"					//  8  Qhi
	  "\t%u"					//  9  Tlo
	  "\t%u"					// 10  Thi
	  ,

	  QLabel,					//  1  Query label
	  TLabel,					//  2  Target label
	  PctId,					//  3  %id
	  AlnLength,				//  4  Alignment length (cols)
	  MismatchCount,			//  5  Mismatches
	  OpenCount,				//  6  Opens
	  QLo,						//  7  Qlo
	  QHi,						//  8  Qhi
	  TLo,						//  9  Tlo
	  THi						// 10  Thi
	  );

	if (AR->IsLocal())
		{
		double BitScore = AR->GetBitScore();
		double Evalue = AR->GetEvalue();
		fprintf(m_fBlast6, "\t%.2g\t%.1f\n",
		  Evalue, BitScore);	// 11  Evalue, 12 BitScore
		}
	else
		fprintf(m_fBlast6, "\t*\t*\n");
	}

void OutputSink::OutputBlast6NoHits(SeqInfo *Query)
	{
	if (m_fBlast6 == 0)
		return;

	const char *QueryLabel = Query->m_Label;

	fprintf(m_fBlast6,
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
	  "\t*"						// 11  Evalue
	  "\t0\n",					// 12  BitScore
	  QueryLabel);
	}
