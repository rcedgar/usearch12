#include "myutils.h"
#include "constaxsink.h"
#include "alignresult.h"
#include "hitmgr.h"
#include "seqinfo.h"
#include "tax.h"
#include "sort.h"
#include <map>

FILE *ConsTaxSink::m_fTab;

ConsTaxSink::ConsTaxSink(bool Local, bool QueryNucleo, bool TargetNucleo)
  : HitSink(Local, QueryNucleo, TargetNucleo)
	{
	m_Maj = (float) oget_flt(OPT_maj);

	LOCK_CLASS();
	if (ofilled(OPT_tabbedout) && m_fTab == 0)
		m_fTab = CreateStdioFile(oget_str(OPT_tabbedout));
	UNLOCK_CLASS();
	}

void ConsTaxSink::OnAllDone()
	{
	CloseStdioFile(m_fTab);
	}

void ConsTaxSink::OnQueryDone(SeqInfo *Query, HitMgr *HM)
	{
	vector<string> TaxStrs;
	const unsigned HitCount = HM->GetHitCount();
	map<string, unsigned> NameToCount;
	vector<string> Names;
	for (unsigned HitIndex = 0; HitIndex < HitCount; ++HitIndex)
		{
		const AlignResult *AR = HM->GetHit(HitIndex);
		const string &TargetLabel = string(AR->m_Target->m_Label);
		string TaxStr;
		GetTaxStrFromLabel(TargetLabel, TaxStr);
		TaxStrs.push_back(TaxStr);

		GetNamesFromTaxStr(TaxStr, Names);
		const unsigned n = SIZE(Names);
		for (unsigned i = 0; i < n; ++i)
			{
			const string &Name = Names[i];
			IncCountMap(NameToCount, Name);
			}
		}

	const unsigned RankCount = GetRankCount();
	bool Found = false;
	string ConsName;
	for (unsigned k = 0; k < RankCount; ++k)
		{
		char Rank = GetRank(RankCount - k - 1);
		for (map<string, unsigned>::const_iterator p = NameToCount.begin();
		  p != NameToCount.end(); ++p)
			{
			const string &Name = p->first;
			if (Name[0] != Rank)
				continue;
			unsigned Count = p->second;
			float Ratio = float(Count)/HitCount;
			if (Ratio >= m_Maj)
				{
				ConsName = Name;
				Found = true;
				break;
				}
			}
		if (Found)
			break;
		}

	string ConsTax = "*";
	if (ConsName != "")
		{
		const unsigned n = SIZE(TaxStrs);
		for (unsigned i = 0; i < n; ++i)
			{
			const string &TaxStr = TaxStrs[i];
			if (NameIsInTaxStr(TaxStr, ConsName))
				{
				bool Ok = TruncateTaxStrAtName(TaxStr, ConsName, ConsTax);
				asserta(Ok);
				break;
				}
			}
		}

	const char *QueryLabel = Query->m_Label;
	LOCK_CLASS();
	fprintf(m_fTab, "%s\t%s\n", QueryLabel, ConsTax.c_str());
	UNLOCK_CLASS();
	}
