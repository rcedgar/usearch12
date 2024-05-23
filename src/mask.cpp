#include "myutils.h"
#include "fastaseqsource.h"
#include "gobuff.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "filetype.h"
#include "mask.h"

MASK_TYPE StrToMaskType(const char *s, MASK_TYPE Default)
	{
	if (s == 0 || s[0] == 0)
		s = "Default";

	if (_stricmp(s, "Default") == 0) return Default;

#define c(x)	if (!_stricmp(s, #x)) return MT_##x;
	c(None)
	c(Default)
	c(Seg)
	c(Dust)
	c(FastNucleo)
	c(FastAmino)
	c(User)
#undef c
	Die("Invalid mask type '%s'", s);
	return MT_Undefined;
	}

const char *MaskTypeToStr(MASK_TYPE Type)
	{
	switch (Type)
		{
#define c(x)	case MT_##x: return #x;
	c(None)
	c(Default)
	c(Seg)
	c(Dust)
	c(FastNucleo)
	c(FastAmino)
	c(User)
#undef c
		}
	Die("Invalid mask type %u", Type);
	return "??";
	}

void MaskSeq(const byte *Seq, unsigned L, MASK_TYPE Type,
  byte *MaskedSeq)
	{
	switch (Type)
		{
	case MT_None:
		{
		for (unsigned i = 0; i < L; ++i)
			MaskedSeq[i] = toupper(Seq[i]);
		return;
		}

	case MT_Seg:
		SegMaskSeq(Seq, L, MaskedSeq);
		return;

	case MT_Dust:
		DustMaskSeq(Seq, L, MaskedSeq);
		return;

	case MT_FastNucleo:
		FastMaskSeq(Seq, L, MaskedSeq, true);
		return;

	case MT_FastAmino:
		FastMaskSeq(Seq, L, MaskedSeq, false);
		return;

	case MT_User:
		if (MaskedSeq != Seq)
			memcpy(MaskedSeq, Seq, L);
		return;

	default:
		Die("MaskSeq(%u)", Type);
		}
	}
