#ifndef syncmer_h
#define syncmer_h

bool IsSyncmer(uint32 KmerCode, uint k, uint s, uint t, uint d,
  bool Open, bool Rotate);
bool IsSyncmer_Seq(const byte *Seq, uint k, uint s, uint t, uint d,
  bool Open, bool Rotate);
uint32 SeqToKmerCode(const byte *Seq, uint k);
const byte *KmerCodeToSeq(uint32 KmerCode, uint k, byte *Seq);

#endif // syncmer_h
