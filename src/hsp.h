#ifndef hsp_h
#define hsp_h	1

struct HSPData
	{
	unsigned Loi;
	unsigned Loj;
	unsigned Leni;
	unsigned Lenj;
	float Score;
	unsigned User;
	
	void Copy(const HSPData &rhs)
		{
		Loi = rhs.Loi;
		Loj = rhs.Loj;
		Leni = rhs.Leni;
		Lenj = rhs.Lenj;
		Score = rhs.Score;
		User = rhs.User;
		}

	unsigned GetLength() const
		{
		if (Leni != Lenj)
			Die("HSP::GetLength(): Leni %u, Lenj %u, Loi %u, Loj %u, Score %.1f",
			  Leni, Lenj, Loi, Loj, Score);

		return Leni;
		}

	unsigned GetHii() const
		{
		assert(Leni > 0);
		return Loi + Leni - 1;
		}

	unsigned GetMidi() const
		{
		return Loi + Leni/2;
		}

	unsigned GetMidj() const
		{
		return Loj + Lenj/2;
		}

	unsigned GetHij() const
		{
		assert(Lenj > 0);
		return Loj + Lenj - 1;
		}

	bool LeftA() const
		{
		return Loi == 0;
		}

	bool LeftB() const
		{
		return Loj == 0;
		}

	bool RightA(unsigned LA) const
		{
		return Loi + Leni == LA;
		}

	bool RightB(unsigned LB) const
		{
		return Loj + Lenj == LB;
		}

	double OverlapFract(const HSPData &HSP) const
		{
		if (Leni == 0 || Lenj == 0)
			return 0.0;

		unsigned MaxLoi = max(Loi, HSP.Loi);
		unsigned MaxLoj = max(Loj, HSP.Loj);
		unsigned MinHii = min(GetHii(), HSP.GetHii());
		unsigned MinHij = min(GetHij(), HSP.GetHij());

		unsigned Ovi = (MinHii < MaxLoi) ? 0 : MinHii - MaxLoi;
		unsigned Ovj = (MinHij < MaxLoj) ? 0 : MinHij - MaxLoj;

		asserta(Ovi <= Leni && Ovj <= Lenj);
		return double(Ovi*Ovj)/double(Leni*Lenj);
		}

/***
	TermGapLeftA = Loi - Loj
    <------------->
                  Loi       Hii
    ---------------===========----------------- A
                ---===========---- B
                Loj         Hij
                                   <------------->
								   TermGapRightA = LA - Hii - (LB - Hij)
***/

	bool IsStaggered(unsigned LA, unsigned LB) const
		{
		int TermGapLeftA = (int) Loi - (int) Loj;
		int TermGapLeftB = (int) Loj - (int) Loi;
		int TermGapRightA = (int) LA - (int) GetHii() - 1 - ((int) LB - (int) GetHij() - 1);
		int TermGapRightB = (int) LB - (int) GetHij() - 1 - ((int) LA - (int) GetHii() - 1);

		if (TermGapLeftA < 0)
			TermGapLeftA = 0;
		if (TermGapLeftB < 0)
			TermGapLeftB = 0;
		if (TermGapRightB < 0)
			TermGapRightB = 0;

		int GapA = TermGapLeftA + TermGapRightA;
		int GapB = TermGapLeftB + TermGapRightB;

		if (GapA == 0 || GapB == 0)
			return false;

		double r = (LA < LB ? double(GapA)/LA : double(GapB)/LB);
		if (r > 0.5)
			return true;
		return false;
		}

	bool operator<(const HSPData &rhs) const
		{
		return Loi < rhs.Loi;
		}

	void LogMe() const
		{
		Log("Loi=%u Loj=%u Li=%u Lj=%u Score=%.1f\n", Loi, Loj, Leni, Lenj, Score);
		}

	void LogMe2() const
		{
		Log("(%u-%u,%u-%u/%.1f)", Loi, GetHii(), Loj, GetHij(), Score);
		}
	};

// Bendpoint
struct BPData
	{
	unsigned Pos;
	bool IsLo;
	unsigned Index;

	void LogMe() const
		{
		Log("BP%s Pos %u Ix %u", (IsLo ? "lo" : "hi"), Pos, Index);
		}
	};

#endif // hsp_h
