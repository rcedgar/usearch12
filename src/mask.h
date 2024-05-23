#ifndef masktype_h
#define masktype_h

#include "myutils.h" // for byte

enum MASK_TYPE
	{
	MT_Undefined,
	MT_Default,
	MT_None,
	MT_Seg,
	MT_Dust,
	MT_FastNucleo,
	MT_FastAmino,
	MT_User,
	};

MASK_TYPE StrToMaskType(const char *s, MASK_TYPE Default);
const char *MaskTypeToStr(MASK_TYPE Type);
void MaskSeq(const byte *Seq, unsigned L, MASK_TYPE Type,
  byte *MaskedSeq);
double GetLowcPct(const byte *Seq, unsigned L);
void SeqToFasta(FILE *f, const byte *Seq, unsigned L, const char *Label = 0);

void SegMask(const byte *Seq, unsigned L, byte *MaskedSeq);
void DustMask(const byte *Seq, unsigned L, byte *MaskedSeq);
void FastMask(const byte *Seq, unsigned L, byte *MaskedSeq);
void SegMaskSeq(const byte *Seq, unsigned L, byte *MaskedSeq);
unsigned DustMaskSeq(const byte *Seq, unsigned L, byte *MaskedSeq);
unsigned FastMaskSeq(const byte *Seq, unsigned L, byte *MaskedSeq, bool Nucleo);

#endif // masktype_h
