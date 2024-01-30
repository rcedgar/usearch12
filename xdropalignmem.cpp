#include "myutils.h"
#include "tracebit.h"
#include "xtype.h"
#include "hsp.h"
#include "alnparams.h"
#include "alnheuristics.h"
#include "xdpmem.h"
#include "objmgr.h"
#include "pathinfo.h"
#include "seqdb.h" // for test only

#define TRACE	0

float XDropFwdFastMem(XDPMem &Mem, const byte *A, unsigned LA, const byte *B, unsigned LB,
  const AlnParams &AP, float X, unsigned &Leni, unsigned &Lenj, PathInfo &PI);

float XDropBwdFastMem(XDPMem &Mem, const byte *A, unsigned LA, const byte *B, unsigned LB,
  const AlnParams &AP, float X, unsigned &Leni, unsigned &Lenj, PathInfo &PI);

float XDropFwdSplit(XDPMem &Mem, const byte *A, unsigned LA, const byte *B, unsigned LB,
  const AlnParams &AP, float X, unsigned &Leni, unsigned &Lenj, PathInfo &PI);

float XDropBwdSplit(XDPMem &Mem, const byte *A, unsigned LA, const byte *B, unsigned LB,
  const AlnParams &AP, float X, unsigned &Leni, unsigned &Lenj, PathInfo &PI);

static float XDropAlignMemMaxL2(XDPMem &Mem, const byte *A, unsigned LA, const byte *B, unsigned LB,
  unsigned AncLoi, unsigned AncLoj, unsigned AncLen, const AlnParams &AP,
  float X, HSPData &HSP, PathInfo &PI, double MaxL2)
	{
#if TRACE
	Log("\n");
	Log("XDropAlignMemMaxL2()\n");
	Log("Anchor A %*.*s %u - %u\n", AncLen, AncLen, A + AncLoi, AncLoi, AncLoi + AncLen - 1);
	Log("Anchor B %*.*s %u - %u\n", AncLen, AncLen, B + AncLoj, AncLoj, AncLoj + AncLen - 1);
	Log(">A;AncLoi=%u;AncLoj=%u;AncLen=%u;\n", AncLoi, AncLoj, AncLen);
	Log("%*.*s\n", LA, LA, A);
	Log(">B\n");
	Log("%*.*s\n", LB, LB, B);
	Log("\n");
#endif
	if (optset_xdrop_save)
		{
		static FILE *f = 0;
		if (f == 0)
			f = CreateStdioFile(opt(xdrop_save));
		fprintf(f, ">A;AncLoi=%u;AncLoj=%u;AncLen=%u;\n", AncLoi, AncLoj, AncLen);
		fprintf(f, "%*.*s\n", LA, LA, A);
		fprintf(f, ">B\n");
		fprintf(f, "%*.*s\n", LB, LB, B);
		fprintf(f, "\n");
		fflush(f);
		}

	if (AncLen <= 1)
		{
		HSP.Score = 0.0f;
		PI.SetEmpty();
		return 0.0f;
		}

	asserta(AncLoi < LA);
	asserta(AncLoj < LB);
	asserta(AncLoi + AncLen <= LA);
	asserta(AncLoj + AncLen <= LB);
	asserta(LA > 0 && LB > 0);

	PI.SetEmpty();

	unsigned AncHii = AncLoi + AncLen - 1;
	unsigned AncHij = AncLoj + AncLen - 1;
	asserta(AncHii < LA);
	asserta(AncHij < LB);

	PathInfo *FwdPI = ObjMgr::GetPathInfo();
	PathInfo *BwdPI = ObjMgr::GetPathInfo();

	const byte *FwdA = A + AncHii;
	const byte *FwdB = B + AncHij;

	unsigned FwdLA = LA - AncHii;
	unsigned FwdLB = LB - AncHij;

	unsigned BwdLeni, BwdLenj;
	float BwdScore = 0.0f;

	if (AncLoi > g_MaxL || AncLoj > g_MaxL)
		{
		BwdScore = XDropBwdSplit(Mem, A, AncLoi+1, B, AncLoj+1, AP, X, BwdLeni, BwdLenj, *BwdPI);
#if	TRACE
		Log("\n");
		Log("XDropBwdSplit:\n");
		LogAln(A + AncLoi + 1, B + AncLoj + 1, BwdPI->GetPath());
#endif
		}
	else
		{
		BwdScore = XDropBwdFastMem(Mem, A, AncLoi+1, B, AncLoj+1, AP, X, BwdLeni, BwdLenj, *BwdPI);
#if	TRACE
		Log("\n");
		Log("XDropBwdFastMem:\n");
		LogAln(A + AncLoi + 1, B + AncLoj + 1, BwdPI->GetPath());
#endif
		}

#if	DEBUG
	{
	unsigned M, D, I;
	BwdPI->GetCounts(M, D, I);
	asserta(M + D == BwdLeni);
	asserta(M + I == BwdLenj);
	}
#endif

	PI.AppendPath(*BwdPI);

	unsigned FwdLeni, FwdLenj;
	float FwdScore = 0.0f;

	if (FwdLA > g_MaxL || FwdLB > g_MaxL)
		{
		FwdScore = XDropFwdSplit(Mem, FwdA, FwdLA, FwdB, FwdLB, AP, X, FwdLeni, FwdLenj, *FwdPI);
#if	TRACE
		Log("\n");
		Log("XDropFwdSplit FwdLA=%u, FwdLB=%u:\n", FwdLA, FwdLB);
		LogAln(FwdA, FwdB, FwdPI->GetPath());
#endif
		}
	else
		{
		FwdScore = XDropFwdFastMem(Mem, FwdA, FwdLA, FwdB, FwdLB, AP, X, FwdLeni, FwdLenj, *FwdPI);
#if	TRACE
		Log("\n");
		Log("XDropFwdFastMem FwdLA=%u, FwdLB=%u:\n", FwdLA, FwdLB);
		Log("A>%*.*s\n", FwdLA, FwdLA, FwdA);
		Log("B>%*.*s\n", FwdLB, FwdLB, FwdB);
		LogAln(FwdA, FwdB, FwdPI->GetPath());
#endif
		}
#if	DEBUG
	{
	unsigned M, D, I;
	FwdPI->GetCounts(M, D, I);
	asserta(M + D == FwdLeni);
	asserta(M + I == FwdLenj);
	}
#endif

// First & last col in anchor are included in xdrop alignmnets
	//for (unsigned i = 2; i < AncLen; ++i)
	//	*PathPtr++ = 'M';

	asserta(AncLen >= 2);
	PI.AppendMs(AncLen-2);

	//for (const char *From = FwdPI->m_Path; *From != 0; ++From)
	//	*PathPtr++ = *From;
	// *PathPtr = 0;

	PI.AppendPath(*FwdPI);

	const float * const *SubstMx = AP.SubstMx;

	float AncScore = 0.0f;
	for (unsigned k = 0; k < AncLen; ++k)
		{
		byte a = A[AncLoi+k];
		byte b = B[AncLoj+k];
		AncScore += SubstMx[a][b];
		}

	byte a = A[AncLoi];
	byte b = B[AncLoj];
	float DupeScore = SubstMx[a][b];
#if	TRACE
	Log("\n");
	Log("Dupe %c%c %.1f\n", a, b, DupeScore);
#endif
	if (AncLen > 1)
		{
		a = A[AncHii];
		b = B[AncHij];
		DupeScore += SubstMx[a][b];
#if	TRACE
		Log("Dupe += %c%c %.1f = %.1f\n", a, b, SubstMx[a][b], DupeScore);
#endif
		}

	HSP.Score = BwdScore + FwdScore + AncScore - DupeScore;
#if	TRACE
	Log("Total = Bwd(%.1f) + Fwd(%.1f) + Anchor(%.1f) - Dupe(%.1f) = %.1f\n",
	  BwdScore, FwdScore, AncScore, DupeScore, HSP.Score);
#endif

	asserta(AncLoi + 1 >= BwdLeni && AncLoj + 1 >= BwdLenj);
	HSP.Loi = AncLoi + 1 - BwdLeni;
	HSP.Loj = AncLoj + 1 - BwdLenj;
	HSP.Leni = BwdLeni + FwdLeni + AncLen - 2;
	HSP.Lenj = BwdLenj + FwdLenj + AncLen - 2;

	ObjMgr::Down(FwdPI);
	ObjMgr::Down(BwdPI);
	FwdPI = 0;
	BwdPI = 0;

#if	DEBUG
	{
	unsigned M, D, I;
	PI.GetCounts(M, D, I);
	asserta(M + D == HSP.Leni);
	asserta(M + I == HSP.Lenj);
	}
#endif

	return HSP.Score;
	}
	
float XDropAlignMem(XDPMem &Mem, const byte *A, unsigned LA, const byte *B, unsigned LB,
  unsigned AncLoi, unsigned AncLoj, unsigned AncLen, const AlnParams &AP,
  float X, HSPData &HSP, PathInfo &PI)
	{
	float ScoreSplit = XDropAlignMemMaxL2(Mem, A, LA, B, LB, AncLoi, AncLoj, AncLen, AP,
	  X, HSP, PI, g_MaxL2);

#if	DEBUG
	{
	unsigned M, D, I;
	PI.GetCounts(M, D, I);
	if (ScoreSplit > 0.0f && (M + D != HSP.Leni || M + I != HSP.Lenj))
		{
		Log("A %5u %*.*s\n", LA, LA, LA, A);
		Log("B %5u %*.*s\n", LB, LB, LB, B);
		Die("PI!=HSP Score %.1f M %u + D %u != HSP.Leni %u || M %u + I %u != HSP.Lenj %u",
		  ScoreSplit, M, D, HSP.Leni, M, I, HSP.Lenj);
		}
	}
#endif
#if	TRACE
	Log("\n");
	Log("XDropAlignMem done, Score = %.1f\n", ScoreSplit);
	LogAln(A + HSP.Loi, B + HSP.Loj, PI.GetPath());
#endif
	return ScoreSplit;
	}

//void Test1()
//	{
//	const char *A = "NIEQFESALKLQSYAGKCPDIEGSGGKTHSKGITHVRNKFLSNAVYESAISLVMNKNKEFYDIFNREIGKKKSNVQAYIAVSKRLLFHIYSIMKNHKPYKEKMVGNYKRYA";
//	const char *B = "NYRDFQTPEQLAKWCGLNPGENESAGKKRKCGITKRGSKYIRVVLVEAAQTIANMKNTGLSRFYKRLSKKKEHNVAIVAVARKLICLIYHLLINQEFYQEVDCRKRKKGRNESCHEPSLKDEHLTDKVAAIVDAFYGMSDSSRKKALLRALEDISVNKPDQKRSSDGGGMSQRLCFKIPRETNTYADVLMALGTASLLNEIYGDEARTTIEDRGDNLSIEVSSDKGLEFLDTFTRLDIGYPFVKQKNDEIIPNGVDDVFDYEKNKEIEEAYNKFIKTSGKRKNKISNKMIDSGFDKPYEPDPNLKLHKILASMRKGWKSDKEFHEYFINNREKISKLAMINLRHLSDPCFEKCSSDDLDKIVSGSQVLLPIGGKGVNRSKPDSTGKSGLPSDFIDWFSEWMKYRGMFKFLLPYRNGDDFKFLVIVPKDISYNSVVSIHNDLFKMNLWGGIKLDIQATLSLAKILVNHSEEYDKEKGSFKMLNRRPNQVIEGLYQAYFKSMGTAASLMNYSFIGLPGWFTIDDNKTAHEFMDILDEHGKCIGALQEEISSDIPLLQDYRNFLSSGDYRHFLEFLSLFGTYVMQRREKNKWVMQFTVQNLRRIFMVDKDYSEIISDEGFLNLASAIRRATVNAQYRKAQGRRDWDIKYGLAQDWKRVAEQPDKLVIAISEFVQQYNAENARHAEENKERRKNITTRDLDQVLNLIKNYGSFLVGMLLLAYGYARDVKETDESNSTEIKEGEIKMLKPSLIDGYASFFKEITSYCPYPYQEKVFSLLDRGKSVVLRAPTGSGKTLSVIMPYLYLRKNNEKMADRIIYALPMRALAFDLYNSTIQTAEKAGFCVVEDPEKRSGEDRDSEICITIQTGELQNDHFFEGDIIFTTIDQLLSGYLTLPVSLPDRLGNINIGAIIGSLIVFDEVHLMEFQRSFLTATEMMNRFEGLAQFFVMTATLTEPTQKWLAKQIKAEPITCSINELKSVPGYEKRSRVYSFENSPISSDTVIKNHRQKTLVICNRVDKAQQIYSELKEKLINADTKIALLHSRFFKRDRSEIEKNLQLWFGKRSEDVNAILVSTQVVEAGIDISCNTMHTELAPANSLIQRAGRCARYSGTGNVFIHEIEDPSNPLPYPKNEVDATRKEIKNFREKELDPDFEKDLVERVHAPLDDSRLIDSLQNRRGEVNEAMDHGHSSSIQKLIREVDSVNIILTPNPHSENLEISWPEMLSVPRSTLYKLKKKGADDWIFKIPIFSDDEEFKGVKWEKVESVRDISWLVALNPKFARYSQEVGLIFGEEGTRCIILNKKSPYKSYKYSCESFQEHTHRVIAEARSLIPLSDKGIMKLSKKMRLDRNHSIEKMILLACALHDTGKLSVRWQGAMNRWQGISNPNAEQFIKGLPLAHTTYDPVKDWIRMQELKIQKGPHAAEGAYAVFTMIFKGSQTMLGDGDDTQNMAIGIMAAIARHHGPRTSNLQPFELIPQAKNVLNDAMKTLEIDLGDMEISETPNKDVIDSFRQYLKISEIDNGIFLYWLSVRIIRLADQRATSKIGAMSRIILSDRGTFLGKSSEQFRITRKDLPDVLVPARGVEQILVLGSGISVSSDAIQMADELGVEVVFASYYGKPMARLIPASLGGTVKTRREQYNAYNDSRGTNLARSFVRGKLKNQASLLKSFAKKWKGERQDLWGEFRNNSNQIEELIPRLNLIKGKVDDVRGEIMGVEGIGAEIYWRTWAKLIPDDWSFPGRDYPAARDQINSLLNFGYYLLEQEVWAAILYAGLDPYAGFLHADRPGKEKLVYDMMEEFRPLVVDRVVVSLAREMRQGHFQDDCRMTKDGIKIASKAFYGRLDERITYREKSQILRNIIRSQASSVATFLRGERPCYEPFTPRWMDTIIIYDISDNSLRARVAKVLLDFGCIRIQKSAFWGIMNHNTREKLRLRLERMMHEKEGNIQFYPMCSKCFSLKDSIGEIYEVEEEDVSVFMKQYLYCPKIIYFDHVLHIPKPPDQKLQTGIEMHDSITAKEKRRKGAIFYDPELDQA";
//	unsigned LA = 111; // (unsigned) strlen(A);
//	unsigned LB = 169; // (unsigned) strlen(B);
//
//	XDPMem Mem;
//	ObjMgr OM;
//	Mem.m_ObjMgr = &OM;
//	PathInfo *PI = Mem.ObjMgr::GetPathInfo();
//	AlnParams AP;
//	AP.InitFromCmdLine(false);
//	AlnHeuristics AH;
//	AH.InitFromCmdLine(AP);
//	float X = AH.XDropG;
//	unsigned Leni, Lenj;
//	float Score = XDropFwdSplit(Mem, (const byte *) A, LA, (const byte *) B, LB, AP, X, Leni, Lenj, *PI);
//
//	Log("A %2u %*.*s\n", LA, LA, LA, A);
//	Log("B %2u %*.*s\n", LB, LB, LB, B);
//	Log("XDropFwdFastMem %.1f, Leni %u, Lenj %u\n", Score, Leni, Lenj);
//	if (Score > 0.0f)
//		LogAln((const byte *) A, (const byte *) B, PI->m_Path);
//	}
//
////static void Test2()
//	{
//	const char *A = "DDLEAFAKQFKQRRIT";
//	const char *B = "DDLEQFAKMFKQRRIV";
//	unsigned LA = (unsigned) strlen(A);
//	unsigned LB = (unsigned) strlen(B);
//
//	XDPMem Mem;
//	ObjMgr OM;
//	Mem.m_ObjMgr = &OM;
//	PathInfo *PI = Mem.ObjMgr::GetPathInfo();
//	AlnParams AP;
//	AP.InitFromCmdLine(false);
//	AlnHeuristics AH;
//	AH.InitFromCmdLine(AP);
//	float X = AH.XDropG;
//	unsigned Leni, Lenj;
//	float Score = XDropBwdFastMem(Mem, (const byte *) A, LA, (const byte *) B, LB, AP, X, Leni, Lenj, *PI);
//
//	Log("A %2u %*.*s\n", LA, LA, LA, A);
//	Log("B %2u %*.*s\n", LB, LB, LB, B);
//	Log("XDropBwdFastMem %.1f, Leni %u, Lenj %u\n", Score, Leni, Lenj);
//	if (Score > 0.0f)
//		LogAln((const byte *) A + LA - Leni, (const byte *) B + LB - Lenj, PI->m_Path);
//	}

void TestXDrop2(const byte *A, unsigned LA, const byte *B, unsigned LB,
  unsigned AncLoi, unsigned AncLoj, unsigned AncLen)
	{
	XDPMem Mem;
	HSPData HSP;
	PathInfo *PI = ObjMgr::GetPathInfo();
	AlnParams AP;
	AP.InitFromCmdLine(false);
	AlnHeuristics AH;
	AH.InitFromCmdLine(AP);
	float X = AH.XDropG;
	Log("A=%*.*s\n", LA, LA, A);
	Log("  ");
	for (unsigned i = 0; i < AncLoi; ++i)
		Log(" ");
	for (unsigned i = 0; i < AncLen; ++i)
		Log("^");
	Log("\n");

	Log("B=%*.*s\n", LB, LB, B);
	Log("  ");
	for (unsigned i = 0; i < AncLoj; ++i)
		Log(" ");
	for (unsigned i = 0; i < AncLen; ++i)
		Log("^");
	Log("\n");

	XDropAlignMem(Mem, (const byte *) A, LA, (const byte *) B, LB, AncLoi, AncLoj, AncLen, AP, X, HSP, *PI);
	Log("A %2u %*.*s\n", LA, LA, LA, A);
	Log("B %2u %*.*s\n", LB, LB, LB, B);
	if (HSP.Score > 0.0f)
		{
		void LogAlnPretty(const byte *A, const byte *B, const char *Path,
		  bool StripTermGaps);
		LogAlnPretty((const byte *) A + HSP.Loi, (const byte *) B + HSP.Loj, PI->GetPath(), true);
		}
	}

unsigned GetIntField(const string &Label, const string &Name)
	{
	vector<string> Fields;
	Split(Label, Fields, ';');
	unsigned N = SIZE(Fields);
	unsigned n = SIZE(Name);
	for (unsigned i = 0; i < N; ++i)
		{
		const string &Field = Fields[i];
		if (Field.substr(0, n) == Name)
			return (unsigned) atoi(Field.substr(n+1, string::npos).c_str());
		}
	Die("Field '%s' not found in label '%s'", Name.c_str(), Label.c_str());
	return 0;
	}

#if	0
void cmd_testxdrop()
	{
	SeqDB Input;
//	Input.FromFasta(opt(testxdrop));
	const unsigned SeqCount = Input.GetSeqCount();
	asserta(SeqCount%2 == 0);

	for (unsigned i = 0; i < SeqCount; i += 2)
		{
		const char *LabelA = Input.GetLabel(i);

		const byte *A = Input.GetSeq(i);
		const byte *B = Input.GetSeq(i+1);

		unsigned LA = Input.GetSeqLength(i);
		unsigned LB = Input.GetSeqLength(i+1);

		unsigned AncLoi = GetIntField(LabelA, "AncLoi");
		unsigned AncLoj = GetIntField(LabelA, "AncLoj");
		unsigned AncLen = GetIntField(LabelA, "AncLen");

		TestXDrop2(A, LA, B, LB, AncLoi, AncLoj, AncLen);
		}
	}
#endif // 0
