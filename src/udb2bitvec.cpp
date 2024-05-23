#include "myutils.h"
#include "udbdata.h"
#include "bitvec.h"

void cmd_udb2bitvec()
	{
	const string &InputFileName = oget_str(OPT_udb2bitvec);
	const string &OutputFileName = oget_str(OPT_output);
	if (InputFileName == "")
		Die("Missing input filename");
	if (OutputFileName == "")
		Die("Missing -output");

	UDBData *udb = new UDBData;
	udb->FromUDBFile(InputFileName);

	const UDBParams &Params = udb->m_Params;
	unsigned SlotCount = Params.m_SlotCount;
	uint32 WordLength = Params.m_WordOnes;
	asserta(Params.m_WordWidth == WordLength);

	BitVec BV;
	BV.Alloc(SlotCount);

	const uint32 *Sizes = udb->m_Sizes;
	unsigned TargetSeqCount = udb->GetSeqCount();
	for (unsigned Word = 0; Word < SlotCount; ++Word)
		{
		unsigned Size = Sizes[Word];
		bool Present = (Size > 0);
		if (Present)
			BV.SetBit(Word);
		}

	for (unsigned Word = 0; Word < SlotCount; ++Word)
		{
		unsigned Size = Sizes[Word];
		bool Present = (Size > 0);
		asserta(Present == BV.GetBit(Word));
		}

	FILE *f = CreateStdioFile(OutputFileName);
	unsigned Bytes = SlotCount/8 + 1;
	WriteStdioFile(f, &WordLength, sizeof(uint32));
	WriteStdioFile(f, BV.m_Vec, Bytes);
	CloseStdioFile(f);
	}
