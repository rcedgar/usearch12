#include "myutils.h"
#include "hspfinder.h"
#include "diagbox.h"
#include "alnheuristics.h"
#include <algorithm>

#define	TRACE	0

unsigned HSPFinder::GetGlobalHSPs(unsigned MinLength, float /* MinFractId */,
  bool StaggerOk, float &HSPFractId)
	{
	const byte *A = m_SA->m_Seq;
	const byte *B = m_SB->m_Seq;

	const unsigned LA = m_SA->m_L;
	const unsigned LB = m_SB->m_L;

	float X = m_AH->XDropGlobalHSP;
	UngappedBlast(X, StaggerOk, MinLength, m_AH->MinGlobalHSPScore);
	Chain();

	unsigned TotalLength = 0;
	unsigned TotalSameCount = 0;

	unsigned N = m_ChainedHSPCount;
#if TRACE
	Log("\n");
	Log("After chain, %u HSPs:\n", m_ChainedHSPCount);
#endif
	for (unsigned i = 0; i < m_ChainedHSPCount; ++i)
		{
		const HSPData &HSP = *m_ChainedHSPs[i];
#if TRACE
		{
		HSP.LogMe();
		for (uint i = HSP.Loi; i <= HSP.GetHii(); ++i)
			Log("%c", A[i]);
		Log("\n");
		for (uint j = HSP.Loj; j <= HSP.GetHij(); ++j)
			Log("%c", B[j]);
		Log("\n");
		}
#endif
		if (HSP.Leni != HSP.Lenj)
			{
			Warning("HSPFinder::GetHSPs, bad HSP");
			HSP.LogMe();

			m_UngappedHSPCount = 0;
			m_ChainedHSPCount = 0;
			return 0;
			}

		TotalLength += HSP.GetLength();
		TotalSameCount += GetHSPIdCount(HSP);
		}

	HSPFractId = TotalLength == 0 ? 0.0f :
	  float(TotalSameCount)/float(TotalLength);
	return m_ChainedHSPCount;
	}
