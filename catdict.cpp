#include "myutils.h"
#include "catdict.h"

unsigned CatDict::GetSampleIndex(const string &SampleName) const
	{
	map<string, unsigned>::const_iterator q =
		m_SampleNameToIndex.find(SampleName);
	if (q == m_SampleNameToIndex.end())
		Die("Sample not found '%s'", SampleName.c_str());
	unsigned SampleIndex = q->second;
	asserta(SampleIndex < SIZE(m_SampleNames));
	return SampleIndex;
	}

unsigned CatDict::GetCatIndex(const string &CatName) const
	{
	map<string, unsigned>::const_iterator r =
		m_CatNameToIndex.find(CatName);
	asserta(r != m_CatNameToIndex.end());
	unsigned CatIndex = r->second;
	asserta(CatIndex < SIZE(m_CatNames));
	return CatIndex;
	}

void CatDict::Init2(OTUTable &OT, const map<string, string> &SampleToCat)
	{
	m_SampleNames = OT.m_SampleNames;
	const unsigned SampleCount = SIZE(m_SampleNames);
	for (unsigned SampleIndex = 0; SampleIndex < SampleCount; ++SampleIndex)
		{
		const string &SampleName = m_SampleNames[SampleIndex];
		if (m_SampleNameToIndex.find(SampleName) !=
		  m_SampleNameToIndex.find(SampleName))
			Die("Duplicate sample name '%s'", SampleName.c_str());
		m_SampleNameToIndex[SampleName] = SampleIndex;
		}

	unsigned CatCount = 0;
	for (map<string, string>::const_iterator p = SampleToCat.begin();
	  p != SampleToCat.end(); ++p)
		{
		const string &SampleName = p->first;
		unsigned SampleIndex = GetSampleIndex(SampleName);

		const string &CatName = p->second;
		if (m_CatNameToIndex.find(CatName) == m_CatNameToIndex.end())
			{
			m_CatNames.push_back(CatName);
			m_CatNameToIndex[CatName] = CatCount;
			++CatCount;
			}
		}

	m_CatIndexToSampleIndexes.resize(CatCount);
	for (map<string, string>::const_iterator p = SampleToCat.begin();
	  p != SampleToCat.end(); ++p)
		{
		const string &SampleName = p->first;
		const string &CatName = p->second;

		unsigned SampleIndex = GetSampleIndex(SampleName);
		unsigned CatIndex = GetCatIndex(CatName);

		m_CatIndexToSampleIndexes[CatIndex].push_back(SampleIndex);
		}
	}

void CatDict::Init1(const map<string, string> &SampleToCat)
	{
	unsigned CatCount = 0;
	unsigned SampleCount = 0;
	for (map<string, string>::const_iterator p = SampleToCat.begin();
	  p != SampleToCat.end(); ++p)
		{
		const string &SampleName = p->first;
		const string &CatName = p->second;

		if (m_SampleNameToIndex.find(SampleName) !=
		  m_SampleNameToIndex.find(SampleName))
			Die("Duplicate sample name '%s'", SampleName.c_str());
		m_SampleNameToIndex[SampleName] = SampleCount;
		m_SampleNames.push_back(SampleName);
		++SampleCount;

		if (m_CatNameToIndex.find(CatName) == m_CatNameToIndex.end())
			{
			m_CatNames.push_back(CatName);
			m_CatNameToIndex[CatName] = CatCount;
			++CatCount;
			}
		}

	m_CatIndexToSampleIndexes.resize(CatCount);
	for (map<string, string>::const_iterator p = SampleToCat.begin();
	  p != SampleToCat.end(); ++p)
		{
		const string &SampleName = p->first;
		const string &CatName = p->second;

		unsigned SampleIndex = GetSampleIndex(SampleName);
		unsigned CatIndex = GetCatIndex(CatName);

		m_CatIndexToSampleIndexes[CatIndex].push_back(SampleIndex);
		}
	}

const char *CatDict::GetCatName(unsigned CatIndex) const
	{
	asserta(CatIndex < SIZE(m_CatNames));
	return m_CatNames[CatIndex].c_str();
	}

const vector<unsigned> &CatDict::GetSampleIndexes(unsigned CatIndex) const
	{
	asserta(CatIndex < SIZE(m_CatIndexToSampleIndexes));
	return m_CatIndexToSampleIndexes[CatIndex];
	}
