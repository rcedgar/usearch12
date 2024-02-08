#include "chunksearcher.h"
#include "globalaligner.h"
#include "hitmgr.h"
#include "outputsink.h"
#include "upclustersink.h"
#include "clustersink.h"
#include "accepter.h"
#include "alnheuristics.h"
#include "alnparams.h"

bool StrandOptToRevComp(bool RequiredOpt, bool Default);

Searcher *MakeClusterSearcher(CMD Cmd, bool Nucleo)
	{
	InitGlobals(Nucleo);

	switch (Cmd)
		{
	case CMD_cluster_otus:
		{
		if (ofilled_flt(OPT_id)) //src_refactor_opts
			Die("-id not supported by cluster_otus");
		oset_fltd(OPT_id, 0);
		break;
		}

	case CMD_cluster_fast:
	case CMD_cluster_smallmem:
		{
		if (!ofilled_flt(OPT_id)) //src_refactor_opts
			Die("Must specify -id");
		break;
		}

	default:
		Die("MakeClusterSearcher(%s)", CmdToStr(Cmd));
		}

	bool AcceptAll = false;
	if (Cmd == CMD_cluster_otus)
		AcceptAll = true;
	Accepter &accepter = *new Accepter(true, AcceptAll);
	Terminator &terminator = *new Terminator(Cmd);

	asserta(!CmdIsLocal(Cmd));

	Aligner *aligner = new GlobalAligner;
	const AlnParams *AP = AlnParams::GetGlobalAP();
	const AlnHeuristics *AH = AlnHeuristics::GetGlobalAH();
	aligner->Init(AP, AH);
	HitMgr &HM = *new HitMgr(0);

	OutputSink &OS = *new OutputSink(false, Nucleo, Nucleo);

	Searcher *searcher = 0;
	SeqDB *seqdb = 0;
	UDBData *udbdata = 0;
	UParseSink *UPS = 0;

// Need UDBUSortedSearcher even if uparse_alignall because ClusterSink needs it.
	UDBUsortedSearcher *US = 0;
	if (Cmd == CMD_cluster_otus)
		{
		if (!ofilled_uns(OPT_maxhits)) //src_refactor_opts
			{
		// hack, not sure what is going on here...
			oset_unsd(OPT_maxhits, 99);
			}
		US = new ChunkSearcher;
		}
	else
		US = new UDBUsortedSearcher;
	US->m_MinFractId = (float) oget_flt(OPT_id); //src_refactor_opts
	UDBParams Params;
	Params.FromCmdLine(Cmd, Nucleo);
	US->m_UDBData->CreateEmpty(Params);
	seqdb = US->m_UDBData->m_SeqDB;

	ObjMgr *OM = ObjMgr::CreateObjMgr();
	US->InitSearcher(&HM, aligner, &accepter, &terminator, OM);
	US->m_RevComp = false;
	udbdata = US->m_UDBData;
	US->m_RevComp = StrandOptToRevComp(false, false);
	searcher = US;

	if (Cmd == CMD_cluster_otus)
		{
		UPClusterSink *CS = new UPClusterSink(Cmd, seqdb, udbdata);
		HM.AddSink(CS);
		}
	else
		{
	// ClusterSink must go before OutputSink, because OS
	// needs cluster index.
		ClusterSink *CS = new ClusterSink(seqdb, udbdata);
		HM.AddSink(CS);
		}

	HM.AddSink(&OS);

	return searcher;
	}
