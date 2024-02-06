#include "myutils.h"
#include "segmask.h"

#define TRACE	0

extern double lnfac[];

const int g_Window = 12;
static int downset, upset;
static const double locut = 2.2;
static const double hicut = 2.5;
static const int maxtrim = 100;
static const bool overlaps = false;
static const int hilenmin = 0;

// One-time init
static int _aaindex[128];
static unsigned char _aaflag[128];
static double *_entray;

// Const'd pointers to one-time init data.
// This to check thread safety, this means there are
// no global variables.
static const int *aaindex = _aaindex;
static const unsigned char *aaflag = _aaflag;
static const double *entray = _entray;

static void stateon(SMSeq *win);
static double entropy(int *sv);
static struct SMSeq *openwin(SMSeq *parent, int start, int length);
static bool shiftwin1(SMSeq *win);
static void closewin(SMSeq *win);

static void appendseg(SMSegment *segs, SMSegment *seg)
  {
  SMSegment *temp;

   temp = segs;
   while (1)
     {
      if (temp->next==NULL)
        {
         temp->next = seg;
         break;
        }
      else
        {
         temp = temp->next;
        }
     }

   return;
  }
static double lnass(int *sv)
	{
	double	ans;
	int	svi, svim1;
	int	Class, total;
	int    i;

	ans = lnfac[20];
	if (sv[0] == 0)
		return ans;

	total = 20;
	Class = 1;
	svim1 = sv[0];
	for (i=0;; svim1 = svi) {
		if (++i==20) {
			ans -= lnfac[Class];
			break;
			}
		else if ((svi = *++sv) == svim1) {
			Class++;
			continue;
			}
		else {
			total -= Class;
			ans -= lnfac[Class];
			if (svi == 0) {
				ans -= lnfac[total];
				break;
				}
			else {
				Class = 1;
				continue;
				}
			}
		}

	return ans;
	}

static double lnperm(int *sv, int tot)
	{
	double ans;
	int i;

	ans = lnfac[tot];
	for (i=0; sv[i]!=0; i++) 
		ans -= lnfac[sv[i]];
	return ans;
	}

static double getprob(int *sv, int total)
  {
	double ans, totseq;

	const double LN20 = 2.9957322735539909;
	totseq = ((double) total) * LN20;

	ans = lnass(sv) + lnperm(sv, total) - totseq;

	return ans;
	}

static void trim(SMSeq *seq, int *leftend, int *rightend)
  {
  struct SMSeq *win;
   double prob, minprob;
   int shift, len, i;
   int lend, rend;
   int minlen;

/* fprintf(stderr, "%d %d\n", *leftend, *rightend);  */
#if	TRACE
   Log("\ntrim()\n");
#endif

   lend = 0;
   rend = seq->length - 1;
   minlen = 1;
   if ((seq->length-maxtrim)>minlen) minlen = seq->length-maxtrim;

   minprob = 1.;
   for (len=seq->length; len>minlen; len--)
     {
      win = openwin(seq, 0, len);
      i = 0;

      shift = true;
      while (shift)
        {
         prob = getprob(win->state, len);
#if	TRACE
		 Log("i=%d prob=%g\n", i, prob);
#endif
         if (prob<minprob)
           {
            minprob = prob;
            lend = i;
            rend = len + i - 1;
#if	TRACE
			Log("  minprob=%g lend=%d rend=%d\n", minprob, lend, rend);
#endif
           }
         shift = shiftwin1(win);
         i++;
#if	TRACE
		 Log("  shift=%d i=%d\n", shift, i);
#endif
        }
      closewin(win);
     }

/* fprintf(stderr, "%d-%d ", *leftend, *rightend);  */

   *leftend = *leftend + lend;
   *rightend = *rightend - (seq->length - rend - 1);

/* fprintf(stderr, "%d-%d\n", *leftend, *rightend);  */

   closewin(seq);
   return;
  }

static int findlo(int i, int limit, double *H)
	{
	int j;

	for (j=i; j>=limit; j--)
		{
		if (H[j]==-1) break;
		if (H[j]>hicut) break;
		}

	return(j+1);
	}

static int findhi(int i, int limit, double *H)
	{
	int j;

	for (j=i; j<=limit; j++)
		{
		if (H[j]==-1) break;
		if (H[j]>hicut) break;
		}

	return(j-1);
	}

static void closewin(SMSeq *win)
	{
	if (win==NULL) return;

	if (win->state!=NULL)       free(win->state);
	if (win->composition!=NULL) free(win->composition);

	free(win);
	}

#ifdef _MSC_VER
static int _cdecl state_cmp(const void *s1, const void *s2)
#else
static int state_cmp(const void *s1, const void *s2)
#endif
	{
	const int *ps1 = (int *) s1;
	const int *ps2 = (int *) s2;
	return *ps2 - *ps1;
	}

static void decrementsv(int *sv, int Class)
	{
	int	svi;

	while ((svi = *sv++) != 0)
		{
		if (svi == Class && *sv < Class)
			{
			sv[-1] = svi - 1;
			break;
			}
		}
	}

static void incrementsv(int *sv, int Class)
	{
	for (;;)
		{
		if (*sv++ == Class)
			{
			sv[-1]++;
			break;
			}
		}
	}

static bool shiftwin1(SMSeq *win)
	{
	int	j, length;
	int	*comp;

	length = win->length;
	comp = win->composition;

	if ((++win->start + length) > win->parent->length)
		{
		--win->start;
#if	TRACE
		Log(" shiftwin1[1] start=%d\n", win->start);
#endif
		return false;
		}

	if (!aaflag[j = win->seq[0]])
		{
		decrementsv(win->state, comp[aaindex[j]]--);
#if	TRACE
		Log(" decsv, j=%d comp[%d]-- = %d\n", j, aaindex[j], comp[aaindex[j]]);
#endif
		}

	j = win->seq[length];
	++win->seq;

	if (!aaflag[j])
		{
		incrementsv(win->state, comp[aaindex[j]]++);
#if TRACE
		Log(" incsv, j=%d comp[%d]++ = %d\n", j, aaindex[j], comp[aaindex[j]]);
#endif
		}

	if (win->entropy > -2.)
		{
		win->entropy = entropy(win->state);
#if	TRACE
		Log(" win->entropy=%g\n", win->entropy);
#endif
		}

	return true;
	}

static bool hasdash(SMSeq *win)
	{
	const char	*seq, *seqmax;
	seq = win->seq;
	seqmax = seq + win->length;
	while (seq < seqmax)
		{
		if (*seq++ == '-')
			return true;
		}
	return false;
	}

static double entropy(int *sv)
	{
	int	*sv0 = sv;
	double	ent;
	int	i, total;
	int	*svmax;
	double	xtotrecip, xsv;

	for (total = 0; (i = *sv) != 0; ++sv)
		total += i;
	svmax = sv;
	ent = 0.0;
	if (total == g_Window)
		{
		for (sv = sv0; sv < svmax; )
			ent += entray[*sv++];
		return ent;
		}
	if (total == 0)
		return 0.0;

	xtotrecip = 1./(double)total;
	for (sv = sv0; sv < svmax; ) {
		xsv = *sv++;
		ent += xsv * log(xsv * xtotrecip);
		}
	return -ent * xtotrecip / LN2;
	}

static bool genwininit()
	{
//	static char _aachar[20];
	const char	*cp, *cp0;
	int		i;
	char	c;

	for (i = 0; i < sizeof(aaindex)/sizeof(aaindex[0]); ++i)
		{
		_aaindex[i] = 20;
		_aaflag[i] = true;
		}

	for (cp = cp0 = "ACDEFGHIKLMNPQRSTVWY"; (c = *cp) != '\0'; ++cp) {
		i = int(cp - cp0);
		_aaindex[c] = i;
		_aaindex[tolower(c)] = i;
	//	_aachar[i] = tolower(c);
		_aaflag[c] = false;
		_aaflag[tolower(c)] = false;
		}

	double	x, xw;

	int Bytes = (g_Window+1) * sizeof(*entray);
	_entray = (double *)malloc(Bytes);
	memset(_entray, 0, Bytes);
	xw = g_Window;
	for (i = 1; i <= g_Window; ++i)
		{
		x = i / xw;
		_entray[i] = -x * log(x) / LN2;
		}

	downset = (g_Window+1)/2 - 1;
	upset = g_Window - downset;

	aaindex = _aaindex;
//	aachar = _aachar;
	aaflag = _aaflag;
	entray = _entray;

	return true;
	}
static bool InitDone = genwininit();

static void enton(SMSeq *win)
	{
	if (win->state==NULL)
		stateon(win);
	win->entropy = entropy(win->state);
	}

static void compon(SMSeq *win)
	{
	int	*comp;
	int	aa;
	const char	*seq, *seqmax;

	win->composition = comp = (int *) calloc(20*sizeof(*comp), 1);
	seq = win->seq;
	seqmax = seq + win->length;

	while (seq < seqmax)
		{
		aa = *seq++;
		if (!aaflag[aa])
			comp[aaindex[aa]]++;
		}
	}

static void stateon(SMSeq *win)
	{
	int	aa, nel, c;

	if (win->composition == NULL)
		compon(win);

	int Bytes = 21*sizeof(win->state[0]);
	win->state = (int *) malloc(Bytes);
	memset(win->state, 0, Bytes);

	for (aa = nel = 0; aa < 20; ++aa)
		{
		if ((c = win->composition[aa]) == 0)
			continue;
		win->state[nel++] = c;
		}
	for (aa = nel; aa < 21; ++aa)
		win->state[aa] = 0;
	assert(nel <= 21);
	qsort(win->state, nel, sizeof(win->state[0]), state_cmp);
	}

static SMSeq *openwin(SMSeq *parent, int start, int length)
	{
	SMSeq *win;

	if (start<0 || length<0 || start+length>parent->length)
		return((SMSeq *) NULL);

	win = (SMSeq *) malloc(sizeof(SMSeq));
	memset(win, 0, sizeof(SMSeq));

	win->parent = parent;
	if (parent->root==NULL)
		{win->root = parent;}
	else
		{win->root = parent->root;}

	win->start = start;
	win->length = length;
	win->seq = parent->seq + start;
	win->punctuation = false;
	win->entropy = -2.0;
	win->state = (int *) NULL;
	win->composition = (int *) NULL;
	stateon(win);
	return win;
	}

static double *seqent(SMSeq *seq)
	{
	SMSeq *win;
	double *H;
	int i, first, last;

	if (g_Window>seq->length)
		return 0;

	H = (double *) malloc((seq->length+1)*sizeof(double));
	for (i=0; i<seq->length; i++)
		H[i] = -1.0;

	win = openwin(seq, 0, g_Window);
	enton(win);

	first = downset;
	last = seq->length - upset;

	for (i=first; i<=last; i++)
		{
		if (seq->punctuation && hasdash(win))
			{
			H[i] = -1;
			shiftwin1(win);
			continue;
			}
		H[i] = win->entropy;
		shiftwin1(win);
		}

	closewin(win);
	return(H);
	}

static void mergesegs(SMSeq *seq, SMSegment *segs)
	{
	SMSegment *seg, *nextseg;
	int len;

	if (overlaps) return;
	if (segs==NULL) return;

	if (segs->begin<hilenmin) segs->begin = 0;

	seg = segs;
	nextseg = seg->next;

	while (nextseg!=NULL)
		{
		if (seg->end>=nextseg->begin)               /* overlapping segments */
			{
			seg->end = nextseg->end;
			seg->next = nextseg->next;
			free(nextseg);
			nextseg = seg->next;
			continue;
			}
		len = nextseg->begin - seg->end - 1;
		if (len<hilenmin)                            /* short hient segment */
			{
			seg->end = nextseg->end;
			seg->next = nextseg->next;
			free(nextseg);
			nextseg = seg->next;
			continue;
			}
		seg = nextseg;
		nextseg = seg->next;
		}

	len = seq->length - seg->end - 1;
	if (len<hilenmin)
		seg->end = seq->length - 1;
	}

static void freesegs(SMSegment *segs)
	{
	SMSegment *temp;
	while (segs!=NULL)
		{
		temp = segs->next;
		free(segs);
		segs = temp;
		}
	}

void segseq(SMSeq *seq, SMSegment **segs, int offset)
	{
	SMSegment *seg, *leftsegs;
	SMSeq *leftseq;
	int first, last, lowlim;
	int loi, hii, i;
	int leftend, rightend, lend, rend;
	double *H;

	H = seqent(seq);
	if (H==NULL)
		return;

#if	TRACE
	Log("\n");
	for (unsigned i = 0; i < seq->length; ++i)
		Log("H[%u] = %g\n", i, H[i]);
#endif

	first = downset;
	last = seq->length - upset;
	lowlim = first;

	for (i=first; i<=last; i++)
		{
		if (H[i]<=locut && H[i]!=-1)
			{
			loi = findlo(i, lowlim, H);
			hii = findhi(i, last, H);

			leftend = loi - downset;
			rightend = hii + upset - 1;

			trim(openwin(seq, leftend, rightend-leftend+1), &leftend, &rightend);

			if (i+upset-1<leftend)   /* check for trigger g_Window in left trim */
				{
				lend = loi - downset;
				rend = leftend - 1;

				leftseq = openwin(seq, lend, rend-lend+1);
				leftsegs = (SMSegment *) NULL;
				segseq(leftseq, &leftsegs, offset+lend);
				if (leftsegs!=NULL)
					{
					if (*segs==NULL)
						*segs = leftsegs;
					else
						appendseg(*segs, leftsegs);
					}
				closewin(leftseq);
				}

			seg = (SMSegment *) malloc(sizeof(SMSegment));
			memset(seg, 0, sizeof(SMSegment));
			seg->begin = leftend + offset;
			seg->end = rightend + offset;
			seg->next = (SMSegment *) NULL;

			if (*segs==NULL)
				*segs = seg;
			else
				appendseg(*segs, seg);

			i = min(hii, rightend+downset);
			lowlim = i + 1;
			}
		}

	free(H);
	mergesegs(seq, *segs);
	}

static void MakeSMSeq(const byte *Seq, unsigned L, SMSeq &S)
	{
	S.start = 0;
	S.parent = 0;
	S.root = 0;
	S.punctuation = false;
	S.state = 0;
	S.entropy = -2.0;
	S.composition = 0;
	S.seq = (const char *) Seq;
	S.length = (int) L;
	}

void SegMaskSeq(const byte *Seq, unsigned L, byte *MaskedSeq)
	{
	SMSeq S;
	MakeSMSeq(Seq, L, S);

	SMSegment *Segs = 0;
	segseq(&S, &Segs, 0);

	for (unsigned i = 0; i < L; ++i)
		MaskedSeq[i] = toupper(Seq[i]);

	for (const SMSegment *Seg = Segs; Seg; Seg = Seg->next)
		{
		unsigned lo = (unsigned) Seg->begin;
		unsigned hi = (unsigned) Seg->end;
		asserta(lo < L);
		asserta(hi < L);
		if (opt(hardmask))
			{
			for (unsigned i = lo; i <= hi; ++i)
				MaskedSeq[i] = 'x';
			}
		else
			{
			for (unsigned i = lo; i <= hi; ++i)
				MaskedSeq[i] = tolower(MaskedSeq[i]);
			}
		}

	freesegs(Segs);
	}
