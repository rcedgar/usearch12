#include "myutils.h"
#include "seqdb.h"
#include "tax.h"
#include "taxy.h"

void cmd_tax_stats()
	{
	SeqDB Input;
	Input.FromFasta(opt(tax_stats));

	Taxy Ty;
	Ty.FromSeqDB(Input);
#if DEBUG
	Ty.Validate();
	Ty.LogMe();
#endif
	}
