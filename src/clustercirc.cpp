#include "myutils.h"
#include "aligncirc.h"
#include "udbusortedsearcher.h"
#include "udbdata.h"
#include "objmgr.h"
#include "hitmgr.h"
#include "clustersink.h"
#include "outputsink.h"
#include "seqsource.h"
#include "filetype.h"
#include "alpha.h"
#include "progress.h"
#include "pathinfo.h"
#include "xdpmem.h"
#include "alignresult.h"

static const unsigned CIRC_BATCH = 128;

unsigned GetMinWindexWordCount(unsigned QueryUniqueWordCount, double FractId,
  unsigned WordLength, bool Nucleo);
bool StrandOptToRevComp(bool RequiredOpt, bool Default);
void InitGlobals(bool Nucleo);

struct CircHit
	{
	bool Hit;
	bool Plus;
	unsigned TargetIndex;
	unsigned WordCount;
	unsigned StartPosQ;
	unsigned TStart;
	unsigned IdCount;
	unsigned ColCount;
	string Path;

	CircHit() { Clear(); }

	void Clear()
		{
		Hit = false;
		Plus = true;
		TargetIndex = UINT_MAX;
		WordCount = 0;
		StartPosQ = 0;
		TStart = 0;
		IdCount = 0;
		ColCount = 0;
		Path.clear();
		}
	};

struct Cand
	{
	unsigned Target;
	unsigned Count;
	bool Plus;
	};

struct WordProbe
	{
	unsigned Top;
	unsigned UniqueCount;
	vector<unsigned> Targets;
	vector<unsigned> Counts;

	WordProbe()
		{
		Top = 0;
		UniqueCount = 0;
		}
	};

struct CircWorker
	{
	UDBUsortedSearcher *US;
	ObjMgr *OM;
	XDPMem Mem;
	PathInfo *PI;
	vector<byte> Q2;

	explicit CircWorker(UDBData *UDB)
		{
		US = new UDBUsortedSearcher(UDB);
		OM = ObjMgr::CreateObjMgr();
		PI = OM->GetPathInfo();
		}
	};

static bool CandBefore(const Cand &a, const Cand &b)
	{
	if (a.Count != b.Count)
		return a.Count > b.Count;
	return a.Target < b.Target;
	}

static bool PassesId(unsigned IdCount, unsigned ColCount, double MinId)
	{
	if (ColCount == 0)
		return false;
	return double(IdCount) / double(ColCount) >= MinId;
	}

static void ClusterCircCB(string &s)
	{
	unsigned C = ClusterSink::GetClusterCount();
	Ps(s, "%s clusters", IntToStr(C));
	}

// Sparse word-to-centroid counts. Only centroids posted under a query word
// are touched. Slots that were incremented are cleared before return.
static void ProbeIndex(UDBUsortedSearcher *US, SeqInfo *Query, WordProbe &P)
	{
	P.Top = 0;
	P.UniqueCount = 0;
	P.Targets.clear();
	P.Counts.clear();

	US->m_Query = Query;
	US->SetQueryWordsAllNoBad();
	US->SetQueryUniqueWords();
	P.UniqueCount = US->m_QueryUniqueWords.Size;

	const unsigned SeqCount = US->GetSeqCount();
	if (SeqCount == 0 || P.UniqueCount == 0)
		return;

	if (US->m_U.MaxSize < SeqCount)
		{
		if (US->m_U.Data != 0 && US->m_U.Size == 0)
			US->m_U.Size = US->m_U.MaxSize;
		US->m_U.Alloc(SeqCount);
		US->m_U.Size = US->m_U.MaxSize;
		}

	unsigned *U = US->m_U.Data;
	const uint32 *Sizes = US->m_UDBData->m_Sizes;
	const uint32 * const *Rows = US->m_UDBData->m_UDBRows;
	const uint32 *Words = US->m_QueryUniqueWords.Data;
	const unsigned NW = US->m_QueryUniqueWords.Size;

	for (unsigned i = 0; i < NW; ++i)
		{
		uint32 Word = Words[i];
		asserta(Word < US->m_UDBData->m_SlotCount);
		const unsigned Size = Sizes[Word];
		if (Size == 0)
			continue;
		const uint32 *Row = Rows[Word];
		for (unsigned j = 0; j < Size; ++j)
			{
			unsigned Target = Row[j];
			asserta(Target < SeqCount);
			if (U[Target] == 0)
				P.Targets.push_back(Target);
			++(U[Target]);
			}
		}

	const unsigned N = SIZE(P.Targets);
	P.Counts.resize(N);
	for (unsigned i = 0; i < N; ++i)
		{
		unsigned Target = P.Targets[i];
		unsigned n = U[Target];
		P.Counts[i] = n;
		if (n > P.Top)
			P.Top = n;
		U[Target] = 0;
		}
	}

static void LoadQueryWords(CircWorker &W, SeqInfo *Query, vector<uint32> &Words)
	{
	W.US->m_Query = Query;
	W.US->SetQueryWordsAllNoBad();
	W.US->SetQueryUniqueWords();
	const unsigned N = W.US->m_QueryUniqueWords.Size;
	Words.resize(N);
	if (N > 0)
		memcpy(Words.data(), W.US->m_QueryUniqueWords.Data, N*sizeof(uint32));
	}

static unsigned CountShared(UDBData *UDB, const vector<uint32> &QWords,
  const byte *T, unsigned TL, bool *Seen)
	{
	const unsigned QN = SIZE(QWords);
	if (QN == 0 || TL == 0 || T == 0)
		return 0;

	for (unsigned i = 0; i < QN; ++i)
		Seen[QWords[i]] = true;

	unsigned n = 0;
	UDBParams &Params = UDB->m_Params;
	const unsigned End = Params.GetLastValidWordPos(TL);
	if (End != UINT_MAX)
		{
		unsigned Step = Params.m_DBStep;
		if (Step == 0)
			Step = 1;
		for (unsigned Pos = 0; Pos <= End; Pos += Step)
			{
			uint32 Word = Params.SeqToWord(T + Pos);
			if (Word == BAD_WORD)
				continue;
			if (Seen[Word])
				{
				++n;
				Seen[Word] = false;
				}
			}
		}

	for (unsigned i = 0; i < QN; ++i)
		Seen[QWords[i]] = false;
	return n;
	}

static bool AlignCand(CircWorker &W, SeqInfo *PlusQ, SeqInfo *MinusQ,
  const Cand &C, double MinId, CircHit &Hit)
	{
	SeqInfo *Orient = C.Plus ? PlusQ : MinusQ;
	asserta(Orient != 0);
	if (Orient->m_L == 0)
		return false;

	SeqInfo *Target = W.OM->GetSeqInfo();
	W.US->GetTargetSeqInfo(C.Target, Target);
	if (Target->m_L == 0)
		{
		Target->Down();
		return false;
		}

	CircAln Aln;
	AlignCircular(W.Mem, *W.PI, W.Q2, Orient->m_Seq, Orient->m_L,
	  Target->m_Seq, Target->m_L, Aln);
	Target->Down();

	if (!PassesId(Aln.IdCount, Aln.ColCount, MinId))
		return false;

	Hit.Hit = true;
	Hit.Plus = C.Plus;
	Hit.TargetIndex = C.Target;
	Hit.WordCount = C.Count;
	Hit.StartPosQ = Aln.StartPosQ % Orient->m_L;
	Hit.TStart = Aln.TStart;
	Hit.IdCount = Aln.IdCount;
	Hit.ColCount = Aln.ColCount;
	Hit.Path.swap(Aln.Path);
	return true;
	}

static void SearchFrozen(CircWorker &W, SeqInfo *Query, bool DoRC,
  double MinId, CircHit &Hit)
	{
	Hit.Clear();

	WordProbe PlusP;
	ProbeIndex(W.US, Query, PlusP);

	SeqInfo *RC = 0;
	WordProbe MinusP;
	if (DoRC)
		{
		RC = W.OM->GetSeqInfo();
		Query->GetRevComp(RC);
		ProbeIndex(W.US, RC, MinusP);
		}

	const bool UsePlus = !(DoRC && MinusP.Top > PlusP.Top);
	WordProbe &P = UsePlus ? PlusP : MinusP;
	if (P.Top == 0)
		{
		if (RC != 0)
			RC->Down();
		return;
		}

	const unsigned WordLen = W.US->m_UDBData->m_Params.m_WordOnes;
	const unsigned MinU = GetMinWindexWordCount(P.UniqueCount, MinId, WordLen, true);

	vector<Cand> Cands;
	const unsigned N = SIZE(P.Targets);
	Cands.reserve(N);
	for (unsigned i = 0; i < N; ++i)
		{
		if (P.Counts[i] < MinU)
			continue;
		Cand C;
		C.Target = P.Targets[i];
		C.Count = P.Counts[i];
		C.Plus = UsePlus;
		Cands.push_back(C);
		}
	sort(Cands.begin(), Cands.end(), CandBefore);

	for (unsigned i = 0; i < SIZE(Cands); ++i)
		{
		if (AlignCand(W, Query, RC, Cands[i], MinId, Hit))
			break;
		}

	if (RC != 0)
		RC->Down();
	}

static void SearchRange(unsigned ThreadIndex, unsigned ThreadCount,
  vector<CircWorker *> *Workers, vector<SeqInfo *> *Batch,
  vector<CircHit> *Hits, bool DoRC, double MinId)
	{
	CircWorker &W = *(*Workers)[ThreadIndex];
	const unsigned N = SIZE(*Batch);
	for (unsigned i = ThreadIndex; i < N; i += ThreadCount)
		SearchFrozen(W, (*Batch)[i], DoRC, MinId, (*Hits)[i]);
	}

static bool SearchNewCentroids(CircWorker &W, SeqInfo *Query, SeqInfo *RC,
  const vector<unsigned> &NewIdx, const CircHit &Frozen, double MinId,
  bool *Seen, CircHit &Hit)
	{
	Hit.Clear();
	if (NewIdx.empty())
		return false;

	vector<uint32> PlusWords;
	vector<uint32> MinusWords;
	LoadQueryWords(W, Query, PlusWords);
	if (RC != 0)
		LoadQueryWords(W, RC, MinusWords);

	UDBData *UDB = W.US->m_UDBData;
	const unsigned WordLen = UDB->m_Params.m_WordOnes;
	const unsigned MinUPlus =
	  GetMinWindexWordCount(SIZE(PlusWords), MinId, WordLen, true);
	const unsigned MinUMinus =
	  GetMinWindexWordCount(SIZE(MinusWords), MinId, WordLen, true);

	vector<Cand> Cands;
	for (unsigned i = 0; i < SIZE(NewIdx); ++i)
		{
		const unsigned Target = NewIdx[i];
		const byte *T = UDB->m_SeqDB->GetSeq(Target);
		const unsigned TL = UDB->m_SeqDB->GetSeqLength(Target);
		const unsigned PlusCount = CountShared(UDB, PlusWords, T, TL, Seen);
		unsigned MinusCount = 0;
		if (RC != 0)
			MinusCount = CountShared(UDB, MinusWords, T, TL, Seen);

		bool Plus = true;
		unsigned Count = PlusCount;
		unsigned MinU = MinUPlus;
		if (MinusCount > PlusCount)
			{
			Plus = false;
			Count = MinusCount;
			MinU = MinUMinus;
			}
		if (Count == 0 || Count < MinU)
			continue;
		if (Frozen.Hit && Count <= Frozen.WordCount)
			continue;

		Cand C;
		C.Target = Target;
		C.Count = Count;
		C.Plus = Plus;
		Cands.push_back(C);
		}

	sort(Cands.begin(), Cands.end(), CandBefore);
	for (unsigned i = 0; i < SIZE(Cands); ++i)
		{
		if (AlignCand(W, Query, RC, Cands[i], MinId, Hit))
			return true;
		}
	return false;
	}

static void AcceptHit(HitMgr *HM, ObjMgr *OM, SeqInfo *Query, SeqInfo *RC,
  UDBData *UDB, const CircHit &Hit)
	{
	asserta(Hit.Hit);
	asserta(Hit.ColCount > 0);
	asserta(SIZE(Hit.Path) == Hit.ColCount);

	SeqInfo *QOrient = Query;
	if (!Hit.Plus)
		{
		asserta(RC != 0);
		QOrient = RC;
		}

	SeqInfo *Target = OM->GetSeqInfo();
	UDB->GetTargetSeqInfo(Hit.TargetIndex, Target);

	PathInfo *PI = OM->GetPathInfo();
	PI->SetEmpty();
	for (unsigned i = 0; i < Hit.ColCount; ++i)
		PI->AppendChar(Hit.Path[i]);

	AlignResult *AR = OM->GetAlignResult();
	AR->CreateGlobal(*QOrient, *Target, *PI, true);
	AR->m_IdCount = Hit.IdCount;
	AR->m_AlnLength = Hit.ColCount;
	AR->m_MismatchCount = 0;
	AR->m_FirstMCol = 0;
	AR->m_LastMCol = Hit.ColCount - 1;
	AR->m_FirstMQPos = 0;
	AR->m_FirstMTPos = 0;
	AR->m_Filled = true;

	const unsigned QL = QOrient->m_L;
	const unsigned Rot = (QL == 0 ? 0 : Hit.StartPosQ % QL);
	AR->m_HSP.Loi = Rot;
	AR->m_HSP.Leni = 1;

	const unsigned TL = Target->m_L;
	unsigned TLo = Hit.TStart;
	if (TL == 0)
		TLo = 0;
	else if (TLo >= TL)
		TLo = TL - 1;
	AR->m_HSP.Loj = TLo;
	AR->m_HSP.Lenj = 1;

	HM->AppendHit(AR);
	AR->Down();
	PI->Down();
	Target->Down();
	}

static void CommitBatch(CircWorker &Main, HitMgr *HM, ObjMgr *OM,
  vector<SeqInfo *> &Batch, vector<CircHit> &Hits, bool DoRC, double MinId,
  bool *Seen)
	{
	UDBData *UDB = Main.US->m_UDBData;
	vector<unsigned> NewIdx;
	NewIdx.reserve(SIZE(Batch));

	for (unsigned i = 0; i < SIZE(Batch); ++i)
		{
		SeqInfo *Query = Batch[i];
		SeqInfo *RC = 0;
		if (DoRC)
			{
			RC = OM->GetSeqInfo();
			Query->GetRevComp(RC);
			}

		CircHit Chosen;
		const bool GotNew = SearchNewCentroids(Main, Query, RC, NewIdx,
		  Hits[i], MinId, Seen, Chosen);

		HM->SetQuery(Query);
		if (GotNew)
			AcceptHit(HM, OM, Query, RC, UDB, Chosen);
		else if (Hits[i].Hit)
			AcceptHit(HM, OM, Query, RC, UDB, Hits[i]);
		HM->OnQueryDone(Query);

		if (!GotNew && !Hits[i].Hit)
			NewIdx.push_back(HM->m_QueryClusterIndex);

		if (RC != 0)
			RC->Down();
		Query->Down();
		}
	}

void cmd_cluster_circ()
	{
	const string &QueryFileName = oget_str(OPT_cluster_circ);
	if (!ofilled(OPT_id))
		Die("Must set -id");
	const double MinId = oget_flt(OPT_id);
	if (MinId < 0.0 || MinId > 1.0)
		Die("-id out of range, should be 0.0 to 1.0");
	if (ofilled(OPT_sort) && oget_str(OPT_sort) != "")
		Die("-cluster_circ keeps input order");
	if (ofilled(OPT_clusters))
		Die("-cluster_circ does not support -clusters");
	if (ofilled(OPT_constax_report))
		Die("-cluster_circ does not support -constax_report");

	bool Nucleo = false;
	GetFileType(QueryFileName, &Nucleo);
	if (!Nucleo)
		Die("-cluster_circ requires nucleotide sequences");

	InitGlobals(true);
	InitCircAligner();

	const bool DoRC = StrandOptToRevComp(false, true);

	SeqSource *SS = MakeSeqSource(QueryFileName);

	UDBParams Params;
	Params.FromCmdLine(CMD_cluster_circ, true);
	UDBData *UDB = new UDBData;
	UDB->CreateEmpty(Params);
	asserta(!UDB->m_Params.DBIsCoded());
	asserta(!UDB->m_Params.DBIsVarCoded());

	ObjMgr *OM = ObjMgr::CreateObjMgr();
	HitMgr *HM = new HitMgr(0);
	ClusterSink *CS = new ClusterSink(UDB->m_SeqDB, UDB);
	OutputSink *OS = new OutputSink(false, true, true);
	HM->AddSink(CS);
	HM->AddSink(OS);

	const unsigned ThreadCount = GetRequestedThreadCount();
	vector<CircWorker *> Workers;
	Workers.reserve(ThreadCount);
	for (unsigned t = 0; t < ThreadCount; ++t)
		Workers.push_back(new CircWorker(UDB));
	CircWorker Main(UDB);

	bool *Seen = myalloc(bool, UDB->m_SlotCount);
	zero_array(Seen, UDB->m_SlotCount);

	bool Any = false;
	ProgressStartSS(*SS, "UCLUST-circ", ClusterCircCB);
	for (;;)
		{
		if (SS->m_EndOfFile)
			break;

		vector<SeqInfo *> Batch;
		Batch.reserve(CIRC_BATCH);
		for (unsigned i = 0; i < CIRC_BATCH; ++i)
			{
			SeqInfo *Query = OM->GetSeqInfo();
			if (!SS->GetNext(Query))
				{
				Query->Down();
				break;
				}
			Batch.push_back(Query);
			}
		if (Batch.empty())
			break;
		Any = true;

		vector<CircHit> Hits(SIZE(Batch));
		unsigned NT = ThreadCount;
		if (NT > SIZE(Batch))
			NT = SIZE(Batch);

		vector<thread *> Ts;
		Ts.reserve(NT);
		for (unsigned t = 0; t < NT; ++t)
			{
			thread *Th = new thread(SearchRange, t, NT, &Workers, &Batch,
			  &Hits, DoRC, MinId);
			Ts.push_back(Th);
			}
		for (unsigned t = 0; t < NT; ++t)
			{
			Ts[t]->join();
			delete Ts[t];
			}

		CommitBatch(Main, HM, OM, Batch, Hits, DoRC, MinId, Seen);
		}
	ProgressDoneSS();

	if (!Any)
		Die("No sequences in input file");

	ClusterSink::OnAllDone(0, 0);
	}
