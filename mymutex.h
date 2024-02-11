#pragma once

#include <mutex>
#include <string>
#include <map>
#include <chrono>

#define USE_MYMUTEX 0
#if !USE_MYMUTEX

#define MUTEX(varname, strname)	mutex varname

#else

#define MUTEX(varname, strname)	mymutex varname(strname)

class mymutex : public std::mutex
	{
public:
	std::string m_name;

public:
	static std::map<std::string, std::chrono::duration<double> > 
	  mymutex::m_name_to_total;
	static std::map<std::string, unsigned> 
	  mymutex::m_name_to_nrcalls;

private:
	static std::mutex g_mylock;

public:
	mymutex(const std::string &name);
	void lock();

public:
	static void report();

private:
	static void init_name(const string &name);
	};
#endif //1or0
