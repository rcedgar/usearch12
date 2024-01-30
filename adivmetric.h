#ifndef adivmetrics_h
#define adivmetrics_h

enum ADIV_METRIC
	{
#define A(x)	ADIV_##x,
#include "adivs.h"
	};

const unsigned ADIV_COUNT = 0 +
#define A(x)	+1
#include "adivs.h"
	;

ADIV_METRIC StrToADivMetric(const string &Name);
const char *ADivMetricToStr(ADIV_METRIC ADiv);
static const char *ADivMetricToStr(unsigned i) { return ADivMetricToStr((ADIV_METRIC) i); }

#endif // adivmetrics_h
