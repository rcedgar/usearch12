#ifndef roccer_h
#define roccer_h

#include <set>

class Roccer
	{
public:
	static void Not(const vector<bool> &In, vector<bool> &Out);
	static void CatsToBools(const vector<string> &TrueCatNames,
	  const string &PosCatName, vector<bool> &IsPosCats);
	static unsigned GetTrueTotal(const vector<bool> &v);
	static unsigned *GetOrder(const vector<float> &Scores);
	static float GetArea2(float TPR1, float FPR1, float TPR2, float FPR2);
	static float GetAUCPoint(unsigned NT, unsigned NF, unsigned NTP, unsigned NFP);
	static float GetAUCPair(float TPR, float FPR);
	static float GetGini(unsigned NT, unsigned NF);
	static float GetGini2(unsigned NTL, unsigned NFL, unsigned NTR, unsigned NFR);
	static float GetAUC(const vector<float> &TPRs, const vector<float> &FPRs);
	static void GetMaxAUC(const vector<float> &TPRs, const vector<float> &FPRs,
	  const vector<float> &Scores, float &MaxAUC, float &MaxAUCScore);
	static void GetXPRs(const vector<bool> &IsPosCats, const vector<float> &PosCatScores,
	  vector<float> &TPRs, vector<float> &FPRs, vector<float> &XPScores);
	static void GetMinGini(const vector<bool> &IsPosCats, const vector<float> &Scores, float &MinGini, float &MinGiniScore);
	static void ToTabbedFile(FILE *f, const vector<float> &TPRs, const vector<float> &FPRs);
	static void ToTabbedFile3(FILE *f, const vector<float> &TPRs, const vector<float> &FPRs,
	  const vector<float> &XPScores);
	};

#endif // roccer_h
