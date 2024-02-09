#include "myutils.h"
#include "label.h"
#include "tax.h"

void StripAllAnnots(string &Label)
	{
	size_t n = Label.find(';');
	if (n == string::npos || n == 0)
		return;
	Label = Label.substr(0, n);
	}

void GetAllAnnots(const string &Label, string &Annots)
	{
	Annots.clear();
	size_t n = Label.find(';');
	if (n == string::npos || n == 0)
		return;
	unsigned L = SIZE(Label);
	for (unsigned i = unsigned(n) + 1; i < L; ++i)
		{
		char c = Label[i];
		Annots.push_back(c);
		}
	}

const char *GetStrField(const string &Label, const string &NameEq,
  string &Value)
	{
	Value.clear();
	vector<string> Fields;
	Split(Label, Fields, ';');
	const unsigned N = SIZE(Fields);
	for (unsigned i = 0; i < N; ++i)
		{
		const string &Field = Fields[i];
		if (StartsWith(Field, NameEq))
			{
			Value = Field.substr(NameEq.size(), string::npos);
			break;
			}
		}
	return Value.c_str();
	}

void StripAnnot(string &Label, const string &NameEq)
	{
	if (Label.find(NameEq) == string::npos)
		return;

	string NewLabel;
	vector<string> Fields;
	Split(Label, Fields, ';');
	const unsigned N = SIZE(Fields);
	for (unsigned i = 0; i < N; ++i)
		{
		const string &Field = Fields[i];
		if (Field.find(NameEq) == 0)
			continue;
		NewLabel += Field + ";";
		}
	if (NewLabel.find('=') == string::npos)
		{
		Label.clear();
		unsigned n = SIZE(NewLabel);
		for (unsigned i = 0; i + 1 < n; ++i)
			Label.push_back(NewLabel[i]);
		}
	else
		Label = NewLabel;
	}

void AppendIntField(string &Label, const string &NameEq, unsigned Value)
	{
	Psasc(Label, "%s%u", NameEq.c_str(), Value);
	}

void AppendStrField(string &Label, const string &NameEqValue)
	{
	Psasc(Label, "%s", NameEqValue.c_str());
	}

void AppendStrField(string &Label, const string &NameEq, const string &Value)
	{
	Psasc(Label, "%s%s", NameEq.c_str(), Value.c_str());
	}

void AppendSize(string &Label, unsigned Size)
	{
	AppendIntField(Label, "size=", Size);
	}

void AppendTaxStr(string &Label, const string &s)
	{
	AppendStrField(Label, "tax=", s);
	}

void StripSize(string &Label)
	{
	StripAnnot(Label, "size=");
	}

void StripTax(string &Label)
	{
	StripAnnot(Label, "tax=");
	StripAnnot(Label, "taxstr=");
	StripAnnot(Label, "sciname=");
	StripAnnot(Label, "taxlev=");
	size_t n = Label.find("\tRoot;");
	if (n != string::npos)
		Label = Label.substr(0, n);
	}

const char *GetTaxStrFromLabel(const string &Label, string &s)
	{
	GetStrField(Label, "tax=", s);
	return s.c_str();
	}

void ReplaceSizeInLabel(string &Label, unsigned NewSize)
	{
	StripSize(Label);
	AppendSize(Label, NewSize);
	}

unsigned GetIntFieldFromLabel(const string &Label, const string &NameEq, unsigned Default)
	{
	const unsigned n = SIZE(NameEq);
	vector<string> Fields;
	Split(Label, Fields, ';');
	const unsigned N = SIZE(Fields);
	for (unsigned i = 0; i < N; ++i)
		{
		const string &Field = Fields[i];
		if (Field.substr(0, n) == NameEq)
			{
			string a = Field.substr(n, string::npos);
			const char *as = a.c_str();
			if (!IsUintStr(as) && Default == UINT_MAX)
				Die("%s not integer >%s", NameEq.c_str(), Label.c_str());
			if (!IsUintStr(as))
				Die("%s invalid value >%s", NameEq.c_str(), Label.c_str());
			unsigned Value = StrToUint(as);
			return Value;
			}
		}
	if (Default == UINT_MAX)
		Die("%s not found in label >%s", NameEq.c_str(), Label.c_str());
	return Default;
	}

unsigned GetSizeFromLabel(const string &Label, unsigned Default)
	{
	unsigned Size = Default;
	const char *p = strstr(Label.c_str(), ";size=");
	if (p != 0)
		Size = (unsigned) atoi(p+6);
	else if (Default == UINT_MAX)
		Die("Missing size= in >%s", Label.c_str());
	return Size;
	}

unsigned GetTaxIdFromLabel(const string &Label, unsigned Default)
	{
	return GetIntFieldFromLabel(Label, "tax=", Default);
	}

void GetAccFromLabel(const string &Label, string &Acc)
	{
	Acc.clear();
	const unsigned N = SIZE(Label);
	for (unsigned i = 0; i < N; ++i)
		{
		char c = Label[i];
		if (c == ' ' || c == '|' || c == ';')
			{
			if (Acc != "gi")
				return;
			}
		Acc += c;
		}
	}

byte GetSplitRankFromLabel(const string &Label)
	{
	string s;
	GetStrField(Label, "split=", s);
	if (SIZE(s) > 1)
		Die("Invalid split in label >%s", Label.c_str());
	return s.c_str()[0];
	}

void GetOTUNameFromLabel(const string &Label, string &OTUName)
	{
	GetStrField(Label, "otu=", OTUName);
	if (!OTUName.empty())
		return;
	
	GetAccFromLabel(Label, OTUName);
	if (OTUName.empty())
		Die("Empty OTU name in label >%s", Label.c_str());
	}

void GetSampleNameFromLabel(const string &Label, string &SampleName)
	{
	SampleName.clear();

	GetStrField(Label, "sample=", SampleName);
	if (!SampleName.empty())
		return;

	GetStrField(Label, "barcodelabel=", SampleName);
	if (!SampleName.empty())
		return;

	if (ofilled(OPT_sample_delim)) //src_refactor_opts
		{
		const string &d = oget_str(OPT_sample_delim); //src_refactor_opts
		size_t n = Label.find(d);
		if (n == string::npos)
			Die("delim '%s' not found in >%s", d.c_str(), Label.c_str());
		SampleName = Label.substr(0, n);
		return;
		}

	unsigned L = SIZE(Label);
	for (unsigned i = 0; i < L; ++i)
		{
		char c = Label[i];
		if (!isalpha(c) && !isdigit(c) && !(c == '_'))
			return;
		SampleName.push_back(c);
		}
	}

void GetTaxNamesFromLabel(const string &Label, vector<string> &Names)
	{
	string TaxStr;
	GetTaxStrFromLabel(Label, TaxStr);
	GetTaxNamesFromTaxStr(TaxStr, Names);
	}
