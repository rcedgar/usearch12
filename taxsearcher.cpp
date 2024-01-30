#include "myutils.h"
#include "accepter.h"
#include "seqinfo.h"
#include "alignresult.h"
#include "aligner.h"
#include "taxsearcher.h"
#include "hitmgr.h"
#include "taxy.h"
#include "randforest.h"
#include "tax.h"

FILE *TaxSearcher::m_fTab = 0;
FILE *TaxSearcher::m_fFeat = 0;
const Taxy *TaxSearcher::m_Taxy = 0;
char TaxSearcher::m_Rank = 0;
char TaxSearcher::m_ParentRank = 0;
const RandForest *TaxSearcher::m_RF = 0;
unsigned TaxSearcher::m_MaxTryCount = 0;
//const map<string, unsigned> *TaxSearcher::m_NameToDepthFirstIndex = 0;

void TaxSearcher::CloseOutputFiles()
	{
	CloseStdioFile(m_fTab);
	CloseStdioFile(m_fFeat);
	}

void TaxSearcher::SetQueryImpl()
	{
	UDBUsortedSearcher::SetQueryImpl();
	ClearTrain();
	}

void TaxSearcher::Init()
	{
	LOCK_CLASS();
	if (m_Rank != 0)
		{
		UNLOCK_CLASS();
		return;
		}

	m_Rank = 'g';
	m_ParentRank = 'f';
	m_MaxTryCount = 32;

	SeqDB *DB = GetSeqDB();
	asserta(DB != 0);
	Taxy *t = new Taxy;
	t->FromSeqDB(*DB);

	m_Taxy = t;
//	m_NameToDepthFirstIndex = &NameToDepthFirstIndex;

	if (optset_tabbedout)
		m_fTab = CreateStdioFile(opt(tabbedout));
	if (optset_featuresout)
		{
		m_fFeat = CreateStdioFile(opt(featuresout));
		WriteFeaturesHdr(m_fFeat);
		}
	if (optset_forestin)
		{
		RandForest *RF = new RandForest;
		RF->FromTabbedFile(opt(forestin));
		m_RF = RF;
		}
	UNLOCK_CLASS();
	}

void TaxSearcher::GetNameAR(AlignResult *AR, char Rank, string &Name) const
	{
	const char *TargetLabel = AR->m_Target->m_Label;
	GetTaxNameFromLabel(TargetLabel, Rank, Name);
	}

void TaxSearcher::GetNameHit(unsigned HitIndex, char Rank, string &Name) const
	{
	AlignResult *AR = m_HitMgr->GetHit(HitIndex);
	return GetNameAR(AR, Rank, Name);
	}

//unsigned TaxSearcher::GetDepthFirstIndex(const string &Name) const
//	{
//	if (Name.empty())
//		return UINT_MAX;
//	map<string, unsigned>::const_iterator p = m_NameToDepthFirstIndex->find(Name);
//	asserta(p != m_NameToDepthFirstIndex->end());
//	unsigned Index = p->second;
//	return Index;
//	}

void TaxSearcher::OnQueryDoneImpl()
	{
	if (m_RF != 0)
		Classify();
	else
		Train();
	UDBUsortedSearcher::OnQueryDoneImpl();
	}

bool TaxSearcher::Align()
	{
	if (m_RF != 0)
		{
		bool Done = UDBUsortedSearcher::Align();
		return Done;
		}

	bool Done = AlignLo();
	if (Done)
		OnDone();
	return Done;
	}

bool TaxSearcher::AlignLo()
	{
	if (opt(train))
		{
		bool Done = Searcher::Align();
		return Done;
		}

	bool Alignable = m_Accepter->AreAlignable(m_Query, m_Target);
	if (!Alignable)
		return ++m_TryCount >= m_MaxTryCount;

	AlignResult *AR = m_Aligner->Align();
	if (AR == 0)
		return ++m_TryCount >= m_MaxTryCount;

	string Name;
	GetNameAR(AR, m_Rank, Name);
	if (Name.empty())
		return ++m_TryCount >= m_MaxTryCount;

	string ParentName;
	GetNameAR(AR, m_ParentRank, ParentName);

	if (m_AR_Top == 0)
		{
		m_Name_Top = Name;
		m_Name_Parent = ParentName;
		ObjMgr::Up(AR);
		m_AR_Top = AR;
		}

	if (m_AR_NN == 0 && Name != m_Name_Top && !Name.empty())
		{
		m_Name_NN = Name;
		ObjMgr::Up(AR);
		m_AR_NN = AR;
		}

	if (m_AR_PN == 0 && ParentName != m_Name_Parent && !ParentName.empty())
		{
		m_Name_PN = ParentName;
		ObjMgr::Up(AR);
		m_AR_PN = AR;
		return true;
		}

	ObjMgr::Down(AR);
	return ++m_TryCount >= m_MaxTryCount;
	}

void TaxSearcher::MakePredStr(const string &TopHitLabel,
  const vector<float> &Probs, string &PredStr) const
	{
	PredStr.clear();

	if (TopHitLabel.empty())
		{
		PredStr = "*";
		return;
		}

	vector<string> Names;
	GetTaxNamesFromLabel(TopHitLabel, Names);
	const unsigned NameCount = SIZE(Names);

	const unsigned CatCount = m_RF->GetCatCount();
	asserta(SIZE(Probs) == CatCount);
	float SumProb = 0.0f;
	vector<string> Preds;
	for (unsigned k = 0; k < NameCount; ++k)
		{
		unsigned NameIndex = NameCount - k - 1;
		const string &Name = Names[NameIndex];
		asserta(SIZE(Name) > 2 && Name[1] == ':');
		string Pred = Name;

		//for (unsigned CatIndex = 0; CatIndex < CatCount; ++CatIndex)
		//	{
		//	const string &CatName = m_RF->GetCatName(CatIndex);
		//	asserta(SIZE(CatName) == 1);
		//	if (CatName[0] == Name[0])
		//		{
		//		float Prob = Probs[CatIndex];
		//		SumProb += Prob;
		//		string ProbStr;
		//		Ps(ProbStr, "(%.4f)", SumProb);
		//		Pred += ProbStr;
		//		break;
		//		}
		//	}

		if (Name[0] == 'g')
			{
			unsigned TCatIndex = m_RF->GetCatIndex("T");
			float TProb = Probs[TCatIndex];
			string s;
			Ps(s, "(%.4f)", TProb);
			Pred += s;
			}

		Preds.push_back(Pred);
		}

	for (unsigned NameIndex = 0; NameIndex < NameCount; ++NameIndex)
		{
		if (NameIndex > 0)
			PredStr += ",";
		PredStr += Preds[NameCount - NameIndex - 1];
		}
	}

void TaxSearcher::Classify()
	{
	Train1(0);

	vector<float> FeatureValues;
	GetFeatureVec(FeatureValues);

	vector<float> Probs;
	m_RF->Classify(FeatureValues, Probs);

	string TopHitLabel;
	AlignResult *AR = m_HitMgr->GetTopHit();
	if (AR != 0)
		TopHitLabel = string(AR->m_Target->m_Label);

	string PredStr;
	MakePredStr(TopHitLabel, Probs, PredStr);

	if (m_fTab != 0)
		{
		LOCK_CLASS();
		fprintf(m_fTab, "%s\t%s\n", m_Query->m_Label, PredStr.c_str());
		UNLOCK_CLASS();
		}
	}

void TaxSearcher::Train1(unsigned StartHit)
	{
	ClearTrain();
	m_StartHit = StartHit;
	const unsigned HitCount = m_HitMgr->GetHitCount();
	for (unsigned HitIndex = StartHit; HitIndex < HitCount; ++HitIndex)
		{
		AlignResult *AR = m_HitMgr->GetHit(HitIndex);

		string Name;
		GetNameAR(AR, m_Rank, Name);

		string ParentName;
		GetNameAR(AR, m_ParentRank, ParentName);

		if (m_AR_Top == 0)
			{
			m_Name_Top = Name;
			m_Name_Parent = ParentName;
			m_AR_Top = AR;
			}

		if (m_AR_NN == 0 && Name != m_Name_Top && !Name.empty())
			{
			m_Name_NN = Name;
			m_AR_NN = AR;
			}

		if (m_AR_PN == 0 && ParentName != m_Name_Parent && !ParentName.empty())
			{
			m_Name_PN = ParentName;
			m_AR_PN = AR;
			break;
			}
		}
	}

void TaxSearcher::Train()
	{
	unsigned HitCount = m_HitMgr->GetHitCount();
	for (unsigned StartHit = 0; StartHit < HitCount; ++StartHit)
		{
		Train1(StartHit);
		WriteTabbed(m_fTab);
		WriteFeatures(m_fFeat);
		}
	}

void TaxSearcher::OnDone()
	{
	if (opt(train))
		return;

	WriteTabbed(m_fTab);
	WriteFeatures(m_fFeat);

#define d(x)	if (m_AR_##x != 0) { ObjMgr::Down(m_AR_##x); m_AR_##x = 0; }
	d(Top)
	d(NN)
	d(PN)
#undef d
	}

void TaxSearcher::WriteTabbed(FILE *f) const
	{
	if (f == 0)
		return;

	LOCK_CLASS();
	UNLOCK_CLASS();
	}

void TaxSearcher::GetFeature(TSFEAT Feat, TSFEAT_VALUE &Val) const
	{
	switch (Feat)
		{
#define F(x, f, y)	\
	case TSF_##x: \
		{ \
		Val.Type = TFT_Float; \
		Val.FloatValue = Get##f##IdAR(m_AR_##y); \
		return; \
		}

#define dF(x, f, y)	\
	case TSF_##x: \
		{ \
		float KiTop = Get##f##IdAR(m_AR_Top); \
		float Ki = Get##f##IdAR(m_AR_##y); \
		Val.Type = TFT_Float; \
		Val.FloatValue = KiTop - Ki; \
		return; \
		}

#define dF(x, f, y)	\
	case TSF_##x: \
		{ \
		float KiTop = Get##f##IdAR(m_AR_Top); \
		float Ki = Get##f##IdAR(m_AR_##y); \
		Val.Type = TFT_Float; \
		Val.FloatValue = KiTop - Ki; \
		return; \
		}

	F(KiTop, Kmer, Top)
	F(FiTop, Fract, Top)

	F(KiNN, Kmer, NN)
	F(FiNN, Fract, NN)

	F(KiPN, Kmer, PN)
	F(FiPN, Fract, PN)

	dF(dKi, Kmer, NN)
	dF(dFi, Kmer, NN)

	dF(dPKi, Kmer, PN)
	dF(dPFi, Fract, PN)

	case TSF_Sib:
		{
		Val.Type = TFT_Int;
		Val.IntValue = GetSiblingCount(m_Name_Top);
		return;
		}

	case TSF_PSib:
		{
		Val.Type = TFT_Int;
		Val.IntValue = GetSiblingCount(m_Name_Parent);
		return;
		}

	//case TSF_TrueName:
	//	{
	//	Val.Type = TFT_Str;
	//	GetTaxNameFromLabel(string(m_Query->m_Label), m_Rank, Val.StringValue);
	//	return;
	//	}

	//case TSF_TopName:
	//	{
	//	Val.Type = TFT_Str;
	//	unsigned HitCount = m_HitMgr->GetHitCount();
	//	if (m_AR_Top == 0)
	//		Val.StringValue = "*";
	//	else
	//		GetTaxNameFromLabel(string(m_ARTop->m_Target->m_Label), m_Rank, Val.StringValue);
	//	return;
	//	}

	//case TSF_DFI:
	//	{
	//	unsigned DFI = GetDepthFirstIndex(m_Name_Top);
	//	return (float) DFI;
	//	}
		}
	asserta(false);
	}

char TaxSearcher::GetTopNameIsTrueName() const
	{
	string TrueName;
	GetTaxNameFromLabel(string(m_Query->m_Label), m_Rank, TrueName);
	if (TrueName.empty())
		return false;

	string TopName;
	if (m_AR_Top != 0)
		GetTaxNameFromLabel(string(m_AR_Top->m_Target->m_Label), m_Rank, TopName);
	if (TopName == TrueName)
		return 'T';
	else
		return 'F';
	}

char TaxSearcher::GetLCR() const
	{
	byte LCR = 0;
	if (m_AR_Top == 0)
		return 'r';
	const string QueryLabel = string(m_Query->m_Label);
	const string TopLabel = string(m_AR_Top->m_Target->m_Label);
	LCR = GetLCRFromLabels(QueryLabel, TopLabel);
	return LCR;
	}

float TaxSearcher::GetFractIdAR(AlignResult *AR) const
	{
	if (AR == 0)
		return 0.0f;
	return (float) AR->GetFractId();
	}

float TaxSearcher::GetKmerIdAR(AlignResult *AR) const
	{
	if (AR == 0)
		return 0.0f;
	return (float) AR->GetKmerId();
	}

unsigned TaxSearcher::GetSiblingCount(const string &Name) const
	{
	if (Name.empty())
		return 0;
	unsigned n = m_Taxy->GetSiblingCount(Name);
	return n;
	}

void TaxSearcher::GetFeatureVec(vector<float> &Values) const
	{
	Values.clear();
#define	T(x)	\
		{ \
		TSFEAT_VALUE Val; \
		GetFeature(TSF_##x, Val); \
		float v; \
		switch (Val.Type) \
			{ \
		case TFT_Float:	v = Val.FloatValue; break; \
		case TFT_Int:	v = float(Val.IntValue);  break; \
		case TFT_Str:	Die("string feature"); \
		default: asserta(false); \
			} \
		Values.push_back(v); \
		}

#include "tsfeats.h"
	}

void TaxSearcher::WriteFeatures(FILE *f) const
	{
	LOCK_CLASS();

	string Acc;
	GetAccFromLabel(m_Query->m_Label, Acc);
	fprintf(f, "%s", Acc.c_str());

	TSFEAT_VALUE Val;

#define	T(x)	\
		{ \
		GetFeature(TSF_##x, Val); \
		switch (Val.Type) \
			{ \
		case TFT_Float:	fprintf(f, "\t%.3g", Val.FloatValue); break; \
		case TFT_Int:	fprintf(f, "\t%u", Val.IntValue);  break; \
		case TFT_Str:	fprintf(f, "\t%s", Val.StringValue.c_str());  break; \
		default: asserta(false); \
			} \
		}

	//char LCR = GetLCR();
	char Cat = GetTopNameIsTrueName();
#include "tsfeats.h"
	fprintf(f, "\t%c", Cat);
	fprintf(f, "\n");

	UNLOCK_CLASS();
	}

void TaxSearcher::WriteFeaturesHdr(FILE *f) const
	{
	if (f == 0)
		return;

	fprintf(f, "Label");
#define	T(x)	fputc('\t', f); fputs(#x, f);
#include "tsfeats.h"
	fprintf(f, "\tLCR=g");
	fprintf(f, "\n");
	}

const char *TaxSearcher::GetFeatureFmt(TSFEAT Feat) const
	{
	switch (Feat)
		{
	case TSF_Sib:
	case TSF_PSib:
//	case TSF_DFI:
		return "%.0f";
		}
	return "%.3g";
	}
