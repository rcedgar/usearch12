#include "myutils.h"
#include "twobit.h"
#include "seqdb.h"
#include "alphabit.h"
#include "alpha.h"
#include "getticks.h"

void Test_TwoBit_Decode(const string &FileName)
	{
	SeqDB Input;
	Input.FromFasta(FileName);
	const uint SeqCount = Input.GetSeqCount();
	vector<byte *> TwoBitSeqs;
	uint MaxL = 0;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		ProgressStep(SeqIndex, SeqCount, "Encoding");
		const byte *Seq = Input.GetSeq(SeqIndex);
		uint L = Input.GetSeqLength(SeqIndex);
		MaxL = max(L, MaxL);
		byte *Seq_TwoBit = myalloc(byte, L);
		TwoBit_Encode(Seq, L, Seq_TwoBit);
		TwoBitSeqs.push_back(Seq_TwoBit);
		}

	byte *DecodeBuffer = myalloc(byte, 2*MaxL);

	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		ProgressStep(SeqIndex, SeqCount, "Validating methods");

		const char *Label = Input.GetLabel(SeqIndex);
		const byte *Seq = Input.GetSeq(SeqIndex);
		uint L = Input.GetSeqLength(SeqIndex);
		const byte *Seq2 = TwoBitSeqs[SeqIndex];

#define ta(Decode_Func)	\
		for (uint i = L; i < L+32; ++i)	\
			DecodeBuffer[i] = 0xff;	\
		Decode_Func(Seq2, L, DecodeBuffer);	\
		for (uint i = 0; i < L; ++i) \
			if (DecodeBuffer[i] != Seq[i]) \
				Die(#Decode_Func"() mismatch DecodeBuffer[%u] = 0x%02x, Seq=0x%02x", \
				  i, DecodeBuffer[i], Seq[i]); \
		for (uint i = L; i < L+32; ++i) \
			if (DecodeBuffer[i] != 0xff) \
				Die(#Decode_Func"() overflow DecodeBuffer[%u] = 0x%02x, L=%u", \
				  i, DecodeBuffer[i], L);

#define tb(Decode_Func)	\
		for (uint i = L; i < L+32; ++i)	\
			DecodeBuffer[i] = 0xff;	\
		Decode_Func(Seq2, L, DecodeBuffer);	\
		for (uint i = 0; i < L; ++i) \
			if (g_LetterToCharNucleo[DecodeBuffer[i]] != Seq[i]) \
				Die(#Decode_Func"() mismatch DecodeBuffer[%u] = 0x%02x, Seq=0x%02x", \
				  i, DecodeBuffer[i], Seq[i]); \
		for (uint i = L; i < L+32; ++i) \
			if (DecodeBuffer[i] != 0xff) \
				Die(#Decode_Func"() overflow DecodeBuffer[%u] = 0x%02x, L=%u", \
				  i, DecodeBuffer[i], L);

#define ta_off(Decode_Func)	\
		for (uint Lo = 0; Lo < 11; ++Lo) \
			{ \
			for (uint i = L; i < L+32; ++i)	\
				DecodeBuffer[i] = 0xff;	\
			Decode_Func(Seq2, Lo, L, DecodeBuffer);	\
			for (uint i = Lo; i < L-Lo; ++i) \
				if (DecodeBuffer[i] != Seq[Lo+i]) \
					Die(#Decode_Func"() mismatch Lo %u DecodeBuffer[%u] = %c, Seq=%c", \
					  Lo, i, DecodeBuffer[i], Seq[Lo+i]); \
			for (uint i = L; i < L+32; ++i) \
				if (DecodeBuffer[i] != 0xff) \
					Die(#Decode_Func"() overflow DecodeBuffer[%u] = 0x%02x, L=%u", \
					  i, DecodeBuffer[i], L); \
			}

#define tb_off(Decode_Func)	\
		for (uint Lo = 0; Lo < 11; ++Lo) \
			{ \
			for (uint i = L; i < L+32; ++i)	\
				DecodeBuffer[i] = 0xff;	\
			Decode_Func(Seq2, Lo, L, DecodeBuffer);	\
			for (uint i = Lo; i < L-Lo; ++i) \
				if (g_LetterToCharNucleo[DecodeBuffer[i]] != Seq[Lo+i]) \
					Die(#Decode_Func"() mismatch Lo %u DecodeBuffer[%u] = 0x%02x, Seq=0x%02x", \
					  Lo, i, DecodeBuffer[i], Seq[Lo+i]); \
			for (uint i = L; i < L+32; ++i) \
				if (DecodeBuffer[i] != 0xff) \
					Die(#Decode_Func"() overflow DecodeBuffer[%u] = 0x%02x, L=%u", \
					  i, DecodeBuffer[i], L); \
			}

		ta(TwoBit_Decode);
		tb(TwoBit_Decode_LetterCodes);
		ta_off(TwoBit_Decode_Offset);
		tb_off(TwoBit_Decode_LetterCodes_Offset);
#undef ta
#undef tb
#undef ta_off
#undef tb_off
		}

	double FirstT = 0;
#define tt(Decode_Func)	\
	{ \
	double MinT = 0; \
	for (uint Iter = 0; Iter < 10; ++Iter) \
		{ \
		TICKS ticks1 = GetClockTicks(); \
		for (uint k = 0; k < 100; ++k) \
			{ \
			for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex) \
				{ \
				uint L = Input.GetSeqLength(SeqIndex); \
				const byte *Seq2 = TwoBitSeqs[SeqIndex]; \
				Decode_Func(Seq2, L, DecodeBuffer); \
				} \
			} \
		TICKS ticks2 = GetClockTicks(); \
		double T = double(ticks2 - ticks1); \
		if (Iter == 0 || T < MinT) \
			MinT = T; \
		} \
	if (FirstT == 0) \
		FirstT = MinT; \
	ProgressLog("%10.0f  %10.10s  %7.3f  " #Decode_Func "\n", MinT, FloatToStr(MinT), MinT/FirstT); \
	}

#define tt_off(Decode_Func)	\
	{ \
	double MinT = 0; \
	for (uint Iter = 0; Iter < 10; ++Iter) \
		{ \
		TICKS ticks1 = GetClockTicks(); \
		for (uint k = 0; k < 100; ++k) \
			{ \
			for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex) \
				{ \
				uint L = Input.GetSeqLength(SeqIndex); \
				const byte *Seq2 = TwoBitSeqs[SeqIndex]; \
				Decode_Func(Seq2, 0, L, DecodeBuffer); \
				} \
			} \
		TICKS ticks2 = GetClockTicks(); \
		double T = double(ticks2 - ticks1); \
		if (Iter == 0 || T < MinT) \
			MinT = T; \
		} \
	if (FirstT == 0) \
		FirstT = MinT; \
	ProgressLog("%10.0f  %10.10s  %7.3f  " #Decode_Func "\n", MinT, FloatToStr(MinT), MinT/FirstT); \
	}

	tt(TwoBit_Decode);
	tt(TwoBit_Decode_LetterCodes);
	tt_off(TwoBit_Decode_Offset);
	tt_off(TwoBit_Decode_LetterCodes_Offset);
	}
