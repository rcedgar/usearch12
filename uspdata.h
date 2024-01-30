#pragma once

/***
Ungapped Segment Pair
LoQ and LoT are always relative to start of plus strand.
When Plus is false, LoQ is start position of the plus strand
segment of Q that aligns to T (if it is rev-comped).
***/
class USPData
	{
public:
	uint LoQ;
	uint LoT;
	uint Length;
	bool Plus;
	union
		{
		void *User1;
		int Score;
		};
	union
		{
		void *User2;
		bool Trimmed;
		};

public:
	uint GetHiQ() const
		{
		asserta(Length > 0);
		return LoQ + Length - 1;
		}

	uint GetHiT() const
		{
		asserta(Length > 0);
		return LoT + Length - 1;
		}
	};
