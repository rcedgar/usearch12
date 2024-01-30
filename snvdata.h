#ifndef snvdata_h
#define snvdata_h

//    0         1  2       3         4
// chr1    110703  C       T       1|1
// chr1    115271  A       C       0/1
// chr1    950806  C       CAG     1|1
// chr1    1769969 CAAAACAAAAACA   CAAAACA,C       2/1

//  1357265 1|1	A and B both have variant.
//   910584 0|1	A matches ref, B has variant.
//   867070 1|0	A has variant, B matches ref.
//   420768 0/1	One of A,B matches ref, other has variant.
//    13147 1|2	A matches first variant, B matches second variant.
//    13067 2|1	A matches second variant, B matches first variant.
//     9638 2/1	A or B matches second variant, the other (B or A) matches the first variant.
//     4578 1/2	Exactly equivalent to 2/1 above.
//        4 1/0	A or B has the variant, the other (B or A) matches the reference.

struct SNVData
	{
	unsigned SeqIndex;
	uint32 Pos;
	string RefStr;
	vector<string> VarStrs;
	byte iA;
	byte iB;
	bool Phased;

	bool Eq(const SNVData *rhs) const
		{
		return rhs->SeqIndex == SeqIndex &&
		    rhs->Pos == Pos &&
			rhs->RefStr == RefStr &&
			rhs->VarStrs == VarStrs &&
			rhs->iA == iA &&
			rhs->iB == iB &&
			rhs->Phased == Phased;
		}

	void GetVarnt(bool TypeA, string &Varnt) const
		{
		Varnt.clear();
		unsigned StrIndex = (TypeA ? iA : iB);
		if (StrIndex == 0)
			{
			Varnt = RefStr;
			return;
			}
		--StrIndex;
		asserta(StrIndex < SIZE(VarStrs));
		Varnt= VarStrs[StrIndex];
		}

	char GetMut(bool TypeA, string &RefCols, string &VarCols) const
		{
		RefCols.clear();
		VarCols.clear();
		unsigned StrIndex = (TypeA ? iA : iB);
		if (StrIndex == 0)
			{
			RefCols = RefStr;
			VarCols = RefStr;
			return 'M';
			}
		--StrIndex;
		asserta(StrIndex < SIZE(VarStrs));
		const string &v = VarStrs[StrIndex];
		if (!optset_readsimv_null)
			asserta(v != RefStr);
		unsigned Lv = SIZE(v);
		unsigned LR = SIZE(RefStr);
		char Mut = 0;
		if (Lv == LR)
			{
			RefCols = RefStr;
			VarCols = v;
			Mut = 'S';
			}
		else if (Lv < LR)
			{
//			asserta(SIZE(v) == 1);
//			asserta(RefStr[0] == v[0]);
			RefCols = RefStr;
			unsigned d = LR - Lv;
			VarCols = v;
			for (unsigned i = 0; i < d; ++i)
				VarCols.push_back('-');
			Mut = 'D';
			}
		else if (Lv > LR)
			{
//			asserta(SIZE(RefStr) == 1);
//			asserta(RefStr[0] == v[0]);
			RefCols = RefStr;
			unsigned d = Lv - LR;
			for (unsigned i = 0; i < d; ++i)
				RefCols.push_back('-');
			VarCols = v;
			Mut = 'I';
			}
		else
			asserta(false);

//		asserta(SIZE(VarCols) == SIZE(RefCols));
		return Mut;
		}

	void WriteMe(FILE *f) const
		{
		if (f == 0)
			return;
		fprintf(f, "Pos %u Ref %s Var ", Pos, RefStr.c_str());
		const unsigned n = SIZE(VarStrs);
		for (unsigned i = 0; i < n; ++i)
			{
			if (i > 0)
				fprintf(f, ",");
			fprintf(f, "%s", VarStrs[i].c_str());
			}
		fprintf(f, " %u%c%u", iA, Phased ? '|' : '/', iB);
		}
	};

#endif // snvdata_h
