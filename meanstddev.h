#ifndef meanstddev_h
#define meanstddev_h

template<class T> class MeanStdDev
	{
public:
	static void Calc(const vector<T> &v, T &Mean, T &StdDev)
		{
		const unsigned N = SIZE(v);
		T Sum = T(0);
		for (unsigned i = 0; i < N; ++i)
			Sum += v[i];
		Mean = Sum/N;

		T Sum2 = T(0);
		for (unsigned i = 0; i < N; ++i)
			{
			T x = v[i];
			T d = (x - Mean);
			Sum2 += d*d;
			}
		StdDev = (T) sqrt(Sum2/N);
		}
	};

#endif // meanstddev_h

