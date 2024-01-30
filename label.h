#ifndef label_h
#define label_h

const char *GetStrField(const string &Label, const string &NameEq,
  string &Value);
void StripAllAnnots(string &Label);
void StripAnnot(string &Label, const string &NameEq);
void AppendIntField(string &Label, const string &NameEq, unsigned Value);
void AppendStrField(string &Label, const string &NameEqValue);
void AppendStrField(string &Label, const string &NameEq, const string &Value);
void AppendSize(string &Label, unsigned Size);
void AppendTaxStr(string &Label, const string &s);
void StripSize(string &Label);
void StripTax(string &Label);
const char *GetTaxStrFromLabel(const string &Label, string &s);
void ReplaceSizeInLabel(string &Label, unsigned NewSize);
unsigned GetIntFieldFromLabel(const string &Label, const string &NameEq, unsigned Default);
unsigned GetSizeFromLabel(const string &Label, unsigned Default);
unsigned GetTaxIdFromLabel(const string &Label, unsigned Default);
void GetAccFromLabel(const string &Label, string &Acc);
byte GetSplitRankFromLabel(const string &Label);
void GetOTUNameFromLabel(const string &Label, string &OTUName);
void GetSampleNameFromLabel(const string &Label, string &SampleName);
void GetAllAnnots(const string &Label, string &Annots);
bool IlluminaLabelPairMatch(const char *Label1, const char *Label2);
void GetTaxNamesFromLabel(const string &Label, vector<string> &Names);
void GetTaxNameFromLabel(const string &Label, char Rank, string &Name);
char GetLowestRankFromLabel(const string &Label);

#endif // label_h
