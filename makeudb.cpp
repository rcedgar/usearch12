#include "myutils.h"
#include "fastaseqsource.h"
#include "seqinfo.h"
#include "objmgr.h"
#include "udbparams.h"
#include "udbdata.h"
#include "udbfile.h"
#include "seqdb.h"
#include <set>

void MaskDB(SeqDB &DB)
	{
	bool Nucleo = DB.GetIsNucleo();
	MASK_TYPE Default = Nucleo ? MT_FastNucleo : MT_FastAmino;
	MASK_TYPE MaskType = StrToMaskType(sopt(dbmask), Default);
	if (MaskType == MT_Default)
		{
		if (DB.GetIsNucleo())
			MaskType = MT_FastNucleo;
		else
			MaskType = MT_FastAmino;
		}

	DB.Mask(MaskType);
	}

static void MakeUDB(CMD Cmd, const string &InputFileName)
	{
	const string &OutputFileName = opt(output);
	if (InputFileName == "" || OutputFileName == "")
		Die("Missing input or output filename");

	if (optset_accel)
		Die("Use --dbaccel, not --accel");

	if (optset_slots)
		{
		asserta(sizeof(uint32 *) >= sizeof(uint32));
		if (double(sizeof(uint32*))*double(opt(slots)) >= double(UINT_MAX))
			{
			double MaxBytes = double(UINT_MAX) + 1;
			double PtrBytes = double(sizeof(uint32 *));
			double MaxSlots = MaxBytes/PtrBytes;
			Die("-slots %u > max %.0f (%s)", opt(slots), MaxSlots, FloatToStr(MaxSlots));
			}
		}

	SeqDB DB;
	bool AllowGaps = false;
	DB.FromFasta(InputFileName, AllowGaps);
	bool DBIsNucleo = DB.GetIsNucleo();

	UDBParams Params;
	Params.FromCmdLine(Cmd, DBIsNucleo);

	UDBData UDB;
	MaskDB(DB);
	UDB.FromSeqDB(DB, Params);
	UDB.ToUDBFile(OutputFileName);
	}

void cmd_makeudb_usearch()
	{
	MakeUDB(CMD_makeudb_usearch, opt(makeudb_usearch));
	}
