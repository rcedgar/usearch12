#include "myutils.h"
#include <chrono>

extern bool g_AbortProgress;

void ProgressThread()
	{
	for (;;)
		{
		if (g_AbortProgress)
			{
			fprintf(stderr, "ProgressThread == Abort\n");
			break;
			}
		this_thread::sleep_for (std::chrono::seconds(1));
		fprintf(stderr, "ProgressThread\n");
		}
	}
