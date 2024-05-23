#ifndef tax_h
#define tax_h

#include <map>
#include <set>

#define RANKS	"rdkpcofgs"

class Tree;
class SeqDB;

// Needs <map>
void GetDictFromTaxStr(const string &TaxStr, map<char, string> &RankToName);

void GetNamesFromTaxStr(const string &TaxStr, vector<string> &Names);
void GetNameFromTaxStr(const string &TaxStr, char Rank, string &Name);
void GetRanksFromTaxStr(const string &TaxStr, string &Ranks);
void GetNameFromNames(const vector<string> &Names, byte Rank, string &Name);

byte GetLCR(const string &TaxStr1, const string &TaxStr2);
byte GetLCRFromLabels(const string &TaxStr1, const string &TaxStr2);

const char *GetRanks();
unsigned GetRankCount();
unsigned GetRankIndex(char Rank);
char GetRank(unsigned Index);
const char *GetRankName(char Rank);

void TaxNameSetFromLabels(const vector<string> &Labels, set<string> &Names);
void TaxNameSetFromTree(const Tree &T, set<string> &Names);
void TaxNameSetFromSeqDB(const SeqDB &DB, set<string> &Names);
void TaxNameSetFromFasta(const string FileName, set<string> &Names);

void TaxPredCutoff(const string &PredWithScores, float Cutoff, string &Pred);

void GetTaxNamesFromTaxStr(const string &TaxStr, vector<string> &Names);
void ParseTaxNameAndScore(const string &NameAndScore, string &Name, float &Score);
void GetLowestRankFromTaxStr(const string &TaxStr, string &Name);
void StripTaxScores(const string &TaxStrWithScores, string &TaxStr);
bool NameIsInTaxStr(const string &TaxStr, const string &Name);
bool TruncateTaxStrAtName(const string &TaxStr, const string &Name,
  string &TruncTax);

#endif // tax_h
