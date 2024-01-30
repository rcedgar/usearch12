// 1. QUERY
SAMOUT_STR(AR->m_Query->m_Label)

// 2. FLAG
SAMOUT_CHAR('\t')
SAMOUT_UINT(AR->GetSAMBits(HitIndex))

// 3. RNAME
SAMOUT_CHAR('\t')
SAMOUT_STR(AR->m_Target->m_Label)

// 4. POS
SAMOUT_CHAR('\t')
SAMOUT_UINT(AR->m_HSP.Loj + 1)

// 5. MAPQ
SAMOUT_CHAR('\t')
SAMOUT_UINT(GetMapQ())

// 6. CIGAR
SAMOUT_CHAR('\t')
SAMOUT_STR(AR->GetCIGAR())

// 7. RNEXT
SAMOUT_STR("\t*")

// 8. PNEXT
SAMOUT_STR("\t0")

// 9. TLEN, not sure why this is always zero (lastz & bowtie).
SAMOUT_STR("\t0")

// 10. SEQ
SAMOUT_CHAR('\t')
SAMOUT_STRN(AR->GetSAMReadSeq(), AR->GetSAMReadSeqLength())

// 11. QUAL
SAMOUT_STR("\t*")

// Tags
// AS:i:78	XN:i:0	XM:i:3	XO:i:1	XG:i:1	NM:i:4	MD:Z:35^T0G4T7T6	YT:Z:UU
SAMOUT_STR("\tAS:i:")
SAMOUT_INT((int) AR->GetSAMScore())

SAMOUT_STR("\tXN:i:")
SAMOUT_UINT(AR->GetTargetSegWildcardCount())

SAMOUT_STR("\tXM:i:")
SAMOUT_UINT(AR->GetMismatchCount())

SAMOUT_STR("\tXO:i:")
SAMOUT_UINT(AR->GetGapOpenCount())

SAMOUT_STR("\tXG:i:")
SAMOUT_UINT(AR->GetGapOpenCount() + AR->GetGapExtCount())

SAMOUT_STR("\tNM:i:")
SAMOUT_UINT(AR->GetDiffCount())

SAMOUT_STR("\tMD:Z:")
SAMOUT_STR(AR->GetMD())

SAMOUT_STR("\tYT:Z:UU")

#undef SAMOUT_CHAR
#undef SAMOUT_STR
#undef SAMOUT_STRN
#undef SAMOUT_UINT
#undef SAMOUT_INT
