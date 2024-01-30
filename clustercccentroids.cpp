#include "myutils.h"
#include "sort.h"
#include "seqdb.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "alignresult.h"
#include "globalaligner.h"

static void ShuffleStrings(vector<string> &v)
	{
	const unsigned N = SIZE(v);
	for (unsigned i = N - 1; i >= 1; --i)
		{
		unsigned j = randu32()%(i + 1);
		
		string &vi = v[i];
		string &vj = v[j];

		v[i] = vj;
		v[j] = vi;
		}
	}

static void DoCC(FILE *fOut, SeqDB &Input, const vector<string> &Labels, uint N)
	{
	asserta(SIZE(Labels) >= N);
	vector<uint> SeqIndexes;
	for (uint i = 0; i < N; ++i)
		{
		uint SeqIndex = Input.GetSeqIndex(Labels[i]);
		SeqIndexes.push_back(SeqIndex);
		}

	vector<vector<float> > PctIdMx(N);
	for (uint i = 0; i < N; ++i)
		PctIdMx[i].resize(N);

	for (uint i = 0; i < N; ++i)
		{
		uint SeqIndexi = SeqIndexes[i];
		SeqInfo *SIQ = ObjMgr::GetSeqInfo();
		Input.GetSI(SeqIndexi, *SIQ);
		PctIdMx[i][i] = 100;
		for (uint j = i+1; j < N; ++j)
			{
			SeqInfo *SIT = ObjMgr::GetSeqInfo();
			AlignResult *AR = ObjMgr::GetAlignResult();
			uint SeqIndexj = SeqIndexes[j];
			Input.GetSI(SeqIndexj, *SIT);
			bool Ok = GlobalAlign_Easy(*SIQ, *SIT, *AR);
			float PctId = 0;
			if (Ok)
				PctId = (float) AR->GetPctId();
			PctIdMx[i][j] = PctId;
			PctIdMx[j][i] = PctId;

			ObjMgr::Down(SIT);
			ObjMgr::Down(AR);
			}
		ObjMgr::Down(SIQ);
		}

	float BestSum = 0;
	uint Besti = 0;
	for (uint i = 0; i < N; ++i)
		{
		float Sum = 0.0;
		for (uint j = 0; j < N; ++j)
			Sum += PctIdMx[i][j];
		if (Sum > BestSum)
			{
			BestSum = Sum;
			Besti = i;
			}
		}
	uint SeqIndex = SeqIndexes[Besti];
	Input.SeqToFasta(fOut, SeqIndex);
	}

void cmd_cluster_cc_centroids()
	{
	const string FastaFileName = string(opt(cluster_cc_centroids));
	const string CCFileName = string(opt(input));
	const string OutputFileName = string(opt(centroids));
	double MaxDist = DBL_MAX;
	if (optset_maxdist)
		MaxDist = opt(maxdist);
	asserta(MaxDist > 0);

	FILE *fOut = CreateStdioFile(OutputFileName);

	SeqDB Input;
	Input.FromFasta(FastaFileName);
	bool Nuc = Input.GetIsNucleo();
	InitGlobals(Nuc);

	string Line;
	vector<string> Fields;
	FILE *fIn = OpenStdioFile(CCFileName);
	vector<string> CCLines;
	while (ReadLineStdioFile(fIn, Line))
		CCLines.push_back(Line);

	const uint CCCount = SIZE(CCLines);
	for (uint CCIndex = 0; CCIndex < CCCount; ++CCIndex)
		{
		ProgressStep(CCIndex, CCCount, "Processing %u ccs", CCCount);
		const string &Line = CCLines[CCIndex];
		Split(Line, Fields, '\t');
		uint CCIndex2 = StrToUint(Fields[0]);
		uint N = StrToUint(Fields[1]);
		asserta(CCIndex2 == CCIndex);
		asserta(SIZE(Fields) == N + 2);

		vector<string> Labels;
		for (uint i = 0; i < N; ++i)
			Labels.push_back(Fields[i+2]);
		if (N > 100)
			{
			ShuffleStrings(Labels);
			N = 100;
			}
		DoCC(fOut, Input, Labels, N);
		}
	CloseStdioFile(fOut);
	}
