#include "myutils.h"
#include "seqinfo.h"
#include "xdpmem.h"
#include "hspfinder.h"
#include "alnheuristics.h"
#include "alignresult.h"
#include "alnparams.h"
#include "objmgr.h"
#include "pathinfo.h"
#include "hsp.h"
#include "cmd.h"
#include "globalaligner.h"

// S. Uliel, A. Fliess, A. Amir, R. Unger, A simple algorithm for detecting circular
// permutations in proteins (1999)
// https://doi.org/10.1093/bioinformatics/15.11.930

bool GlobalAlign_Circle(XDPMem &Mem, const SeqInfo &Query, const SeqInfo &Target,
  const AlnParams &AP, const AlnHeuristics &AH, HSPFinder &HF,
  PathInfo &PI, bool FullDPAlways, bool FailIfNoHSPs)
	{
	asserta(!FullDPAlways);
	const byte *Q = Query.m_Seq;
	const uint QL = Query.m_L;
	const uint QL2 = 2*QL;

	const byte *T = Target.m_Seq;
	const uint TL = Target.m_L;

// Make Query2=tandem duplication of Query
	SeqInfo *Query2 = ObjMgr::GetSeqInfo();
	PathInfo *PI2 = ObjMgr::GetPathInfo();

	Query2->AllocSeq(QL2);
	Query2->SetLabel("Query2");
	const byte *Q2 = Query2->m_Seq;

	memcpy(Query2->m_SeqBuffer, Q, QL);
	memcpy(Query2->m_SeqBuffer + QL, Q, QL);
	Query2->m_L = 2*QL;

	ViterbiFastMem(Mem, Query2->m_Seq, QL2, Target.m_Seq, TL, AP, *PI2);

	LogAln(Q2, T, PI2->GetPath());
	return false;
	}
