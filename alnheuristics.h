#ifndef alnheuristics_h
#define alnheuristics_h

struct AlnParams;

struct AlnHeuristics
	{
	bool FullDPAlways;
	unsigned BandRadius;
	unsigned HSPFinderWordLength;

	float XDropG;			//  GappedBlast default
	float XDropU;			//  UngappedBlast default
	unsigned MinGlobalHSPLength;

// For HSPs in global alignment
	float XDropGlobalHSP;
	float MinGlobalHSPScore;
	float MinGlobalHSPFractId;

	AlnHeuristics();
	void LogMe() const;
	void InitFromCmdLine(const AlnParams &AP);
	void InitGlobalFullDP();
	static const AlnHeuristics *GetGlobalAH();
	};

#endif // alnheuristics_h
