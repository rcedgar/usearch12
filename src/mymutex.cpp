#include "myutils.h"
#include "mymutex.h"

#if USE_MYMUTEX
std::map<std::string, std::chrono::duration<double> > 
  mymutex::m_name_to_total;
std::map<std::string, unsigned> 
  mymutex::m_name_to_nrcalls;
std::mutex mymutex::g_mylock;

void mymutex::init_name(const string &name)
	{
	g_mylock.lock();
	std::map<std::string, std::chrono::duration<double>>::iterator iter =
	  m_name_to_total.find(name);
	if (iter == m_name_to_total.end())
		{
		m_name_to_total[name] = std::chrono::duration<double>(0);
		m_name_to_nrcalls[name] = 0;
		}
	g_mylock.unlock();
	}

mymutex::mymutex(const std::string &name) : mutex()
	{
	m_name = name;
	}

void mymutex::lock()
	{
	init_name(m_name);
    auto t1 = std::chrono::high_resolution_clock::now();
	mutex::lock();
    auto t2 = std::chrono::high_resolution_clock::now();
	g_mylock.lock();
	std::map<std::string, std::chrono::duration<double> >::iterator iter =
	  m_name_to_total.find(m_name);
	if (iter == m_name_to_total.end())
		{
		g_mylock.unlock();
		fprintf(stderr, "\nmymutex::lock(%s)\n", m_name.c_str());
		exit(1);
		}
	iter->second += t2 - t1;
	m_name_to_nrcalls[m_name] += 1;
	g_mylock.unlock();
	}

void mymutex::report()
	{
	FILE *f = fopen("mymutex.txt", "w");
	for (std::map<std::string, std::chrono::duration<double>>::iterator iter =
	  m_name_to_total.begin(); iter != m_name_to_total.end(); ++iter)
		{
		const string &name = iter->first;
		unsigned nrcalls = m_name_to_nrcalls[name];
		double total = iter->second.count();
		fprintf(f, "%10.4g  %10u  %s\n", total, nrcalls, name.c_str());
		}
	fclose(f);
	}
#endif // USE_MYMUTEX
