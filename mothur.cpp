#include "myutils.h"
#include "otutab.h"

/***
Shared file: mothur OTU table.

label     Group  numOtus  OTU_1   OTU_2  OTU_3  OTU_4  OTU_5   OTU_6  OTU_7   OTU_8   OTU_9
usearch    F3D0        9    749      29    613    426     149    366    196      46      26
usearch    F3D1        9     85       9    441    140     115    372    210      74      14
usearch  F3D141        9    535      57    497    378      38    392    370     169       6
usearch  F3D142        9    313      14    312    255      10    327     92      87       0
usearch  F3D143        9    372       2    247    237      19    185    107     109       1
usearch  F3D144        9    607      14    472    382      25    313     48     171       4
usearch  F3D145        9    849      22    719    627      21    542    155     209       8
usearch  F3D146        9    493      16    349    330      43    248     74     120      18
usearch  F3D147        9   2025     127   1720   1417      31   1367    105     864      11
***/

void OTUTable::ToMothurSharedFile(const string &FileName) const
	{
	if (FileName == "")
		return;

	Progress("Writing %s ...", FileName.c_str());
	FILE *f = CreateStdioFile(FileName);

	fprintf(f, "label\tGroup\tnumOtus");
	unsigned OTUCount = GetOTUCount();
	for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
		fprintf(f, "\t%s", GetOTUName(OTUIndex));
	fprintf(f, "\n");

	unsigned SampleCount = GetSampleCount();
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		fprintf(f, "usearch\t%s", GetSampleName(SampleIndex));
		fprintf(f, "\t%u", OTUCount);
		
		for (unsigned OTUIndex = 0; OTUIndex < OTUCount; ++OTUIndex)
			{
			unsigned Count = GetCount(OTUIndex, SampleIndex);
			fprintf(f, "\t%u", Count);
			}
		fprintf(f, "\n");
		}

	CloseStdioFile(f);
	Progress("done.\n");
	}
