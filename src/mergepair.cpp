#include "myutils.h"
#include "merge.h"

bool MergePair(MergeThreadData &TD)
	{
	bool Ok = MergePre(TD.SI1, true);
	if (!Ok)
		return false;

	Ok = MergePre(TD.SI2, false);
	if (!Ok)
		return false;

	TD.SI2->GetRevComp(TD.SI2RC);

	Ok = MergeAlign(TD);
	if (!Ok)
		return false;

	Ok = MergePost(TD);
	if (!Ok)
		return false;

	return true;
	}
