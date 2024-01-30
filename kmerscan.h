#ifndef kmerscan_h
#define kmerscah_h

#include "onkmer.h"

void KmerScan(const byte *Seq2, uint Lo, uint Len, uint k, fn_OnKmer OnKmer, void *UserData);
void KmerScanBoth(const byte *Seq2, uint Lo, uint Len, uint k, fn_OnKmerBoth OnKmer, void *UserData);
void KmerScan_Skip(const byte *Seq, uint Lo, uint Len, uint k, fn_OnKmer_Skip OnKmer, void *UserData);
void KmerScanBoth_Skip(const byte *Seq2, uint Lo, uint Len, uint k, fn_OnKmerBoth_Skip OnKmer, void *UserData);

void SyncmerScan(const byte *Seq2, uint Lo, uint Len, uint k, uint d, fn_OnKmer OnSyncmer, void *UserData);
void SyncmerScan_Skip(const byte *Seq2, uint Lo, uint Len, uint k, uint d, fn_OnKmer_Skip OnSyncmer, void *UserData);
void SyncmerScanBoth(const byte *Seq2, uint Lo, uint Len, uint k, uint d, fn_OnKmerBoth OnSyncmer, void *UserData);
void SyncmerScanBoth_Skip(const byte *Seq2, uint Lo, uint Len, uint k, uint d, fn_OnKmerBoth_Skip OnSyncmer, void *UserData);

void TwoBit_KmerScan(const byte *Seq2, uint Lo, uint Len, uint k, fn_OnKmer OnKmer, void *UserData);
void TwoBit_KmerScanBoth(const byte *Seq2, uint Lo, uint Len, uint k, fn_OnKmerBoth OnKmer, void *UserData);
void TwoBit_KmerScan_Skip(const byte *Seq, uint Lo, uint Len, uint k, fn_OnKmer_Skip OnKmer, void *UserData);
void TwoBit_KmerScanBoth_Skip(const byte *Seq2, uint Lo, uint Len, uint k, fn_OnKmerBoth_Skip OnKmer, void *UserData);

void TwoBit_SyncmerScan(const byte *Seq2, uint Lo, uint Len, uint k, uint d, fn_OnKmer OnSyncmer, void *UserData);
void TwoBit_SyncmerScan_Skip(const byte *Seq2, uint Lo, uint Len, uint k, uint d, fn_OnKmer_Skip OnSyncmer, void *UserData);
void TwoBit_SyncmerScanBoth(const byte *Seq2, uint Lo, uint Len, uint k, uint d, fn_OnKmerBoth OnSyncmer, void *UserData);
void TwoBit_SyncmerScanBoth_Skip(const byte *Seq2, uint Lo, uint Len, uint k, uint d, fn_OnKmerBoth_Skip OnSyncmer, void *UserData);

#endif // kmerscan_h
