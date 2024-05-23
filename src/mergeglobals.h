#ifndef STORECLASS
#error "STORECLASS"
#endif

//STORECLASS omp_lock_t g_GetNextLock;
//STORECLASS omp_lock_t g_MergeOutLock;
//STORECLASS omp_lock_t g_TotalsLock;
//STORECLASS omp_lock_t g_ReportLock;

STORECLASS FILE *g_fFastqOut;
STORECLASS FILE *g_fFastaOut;
STORECLASS FILE *g_fEEOut;
STORECLASS FILE *g_fFqNotmergedFwd;
STORECLASS FILE *g_fFqNotmergedRev;
STORECLASS FILE *g_fFaNotmergedFwd;
STORECLASS FILE *g_fFaNotmergedRev;
STORECLASS FILE *g_fFqOverlapFwd;
STORECLASS FILE *g_fFqOverlapRev;
STORECLASS FILE *g_fFaOverlapFwd;
STORECLASS FILE *g_fFaOverlapRev;
STORECLASS FILE *g_fRep;
STORECLASS FILE *g_fAln;

STORECLASS unsigned g_TailCount1;
STORECLASS unsigned g_TailCount2;
STORECLASS unsigned g_MergedTooShortCount;
STORECLASS unsigned g_MergedTooLongCount;
STORECLASS unsigned g_OvTooShortCount;
STORECLASS unsigned g_TooShortCount1;
STORECLASS unsigned g_TooShortCount2;
STORECLASS unsigned g_MaxEECount;
STORECLASS unsigned g_StaggeredCount;
STORECLASS unsigned g_MinQCount;
STORECLASS unsigned g_InRecCount;
STORECLASS unsigned g_OutRecCount;
STORECLASS unsigned g_ExactOverlapCount;
STORECLASS unsigned g_NonExactOverlapCount;
STORECLASS unsigned g_NotAlignedCount;
STORECLASS unsigned g_MaxDiffsCount;
STORECLASS double g_SumOvLength;
STORECLASS double g_SumMergedLength;
STORECLASS double g_SumEE1;
STORECLASS double g_SumEE2;
STORECLASS double g_SumMergedEE;
STORECLASS vector<unsigned> *g_MergeLengths;

#undef STORECLASS
