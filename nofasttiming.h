#ifndef nofasttiming_h
#define nofasttiming_h

#undef FAST_TIMING
#undef StartFastTimer
#undef EndFastTimer
#define FAST_TIMING	0

#define InitFastTiming()	/* empty */
#define EndFastTiming()		/* empty */
#define FastTimer(x)		/* empty */

#endif // nofasttiming_h
