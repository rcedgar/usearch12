#pragma once

static inline uint GetLoj(uint i, uint dLr)
	{
	asserta(i > 0);
	int Loj = int(i) - int(dLr);
	if (Loj < 1)
		return 1;
	return uint(Loj);
	}

static inline uint GetHij(uint i, uint dLr, uint LT)
	{
	asserta(i > 0);
	int Hij = int(i) + int(dLr);
	if (Hij > int(LT))
		Hij = int(LT);
	return uint(Hij);
	}

static inline uint GetLoj0(uint i, uint dLr)
	{
	int Loj = int(i) - int(dLr);
	if (Loj < 0)
		return 0;
	return uint(Loj);
	}

static inline uint GetHij0(uint i, uint dLr, uint LT)
	{
	int Hij = int(i) + int(dLr);
	if (Hij >= int(LT))
		Hij = int(LT) - 1;
	return uint(Hij);
	}
