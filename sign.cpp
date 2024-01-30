#include "myutils.h"
#include "sign.h"
#include <time.h>

static char *DecryptSig(const SigData *SD)
	{
	static char Str[MAXSIGSTR+1];
	unsigned CheckSum = 0;
	for (unsigned i = 0; i < MAXSIGSTR; ++i)
		{
		byte c = SD->Str[i];
		CheckSum += c;
		}
	if (CheckSum != SD->CheckSum)
		{
		//byte *p = (byte *) SD;
		//for (;;)
		//	*p++ = 0;
		return 0;
		}

	char cKey = (char) SD->Key;
	for (unsigned i = 0; i < MAXSIGSTR; ++i)
		{
		char c = SD->Str[i]^cKey;
		Str[i] = c;
		if (c == 0)
			break;
		}
	return Str;
	}

const char *GetSig()
	{
	//if (Sigs[0].Key == 0xffffffff)
	//	return "";

	srand((unsigned) time(0));
	int ix = rand()%3;
	const SigData *SD = Sigs + ix;
	const char *Sig = DecryptSig(SD);
	if (Sig == 0)
		exit(23);
	return Sig;
	}

static void PokeSig(SigData *SD, const string &Str)
	{
	unsigned L = SIZE(Str);
	if (L > MAXSIGSTR)
		Die("String too long");

	const char *s = Str.c_str();
	srand((unsigned) time(0));
	SD->Key = (unsigned) rand();
	SD->RandomJunk = (unsigned) rand();

	char cKey = (char) SD->Key;
	for (unsigned i = 0; i < L; ++i)
		SD->Str[i] = s[i]^cKey;
	SD->Str[L] = 0^cKey;
	for (unsigned i = L+1; i < MAXSIGSTR; ++i)
		SD->Str[i] = (char) rand();

	unsigned CheckSum = 0;
	for (unsigned i = 0; i < MAXSIGSTR; ++i)
		{
		byte c = SD->Str[i];
		CheckSum += c;
		}
	SD->CheckSum = CheckSum;

	if (opt(verbose))
		{
		unsigned n = (unsigned) sizeof(SigData);
		printf("byte SigDataBuffer[%u] =\n", n);
		printf("\t{\n");
		const byte *Buff = (const byte *) SD;
		for (unsigned i = 0; i < n; ++i)
			{
			if (i > 0 && i%10 == 0)
				printf("\n\t");
			printf("0x%02x,", Buff[i]);
			}
		printf("\t};\n");
		}
	}

void Sign(const string &InputFileName, const string &OutputFileName,
  const string &Str)
	{
	if (InputFileName == "")
		Die("Missing --input");
	if (OutputFileName != "" && Str == "")
		Die("Missing --str");

	uint64 otFileSize;
	byte *Buffer = ReadAllStdioFile64(InputFileName, otFileSize);
	unsigned FileSize = (unsigned) otFileSize;
	if ((off_t) FileSize != otFileSize)
		Die("File too big");

	byte *Start = Buffer;
	byte *End = Buffer + FileSize - sizeof(SigData);
	unsigned Count = 0;
	Progress("OutputFileName='%s'\n", OutputFileName.c_str());
	for (byte *p = Start; p < End; ++p)
		{
		SigData *SD = (SigData *) p;
		if (SD->Magic1 == MAGIC1 && SD->Magic2 == MAGIC2)
			{
			unsigned Offset = (unsigned) (p - Buffer);
			if (OutputFileName.empty())
				{
				const char *StrN = DecryptSig(SD);
				Progress("Sig #%u found at %u, str=\"%s\"\n", ++Count, Offset, StrN);
				}
			else
				{
				Progress("Write sig #%u found at %u\n", ++Count, Offset);
				PokeSig(SD, Str);
				}
			}
		}
	Progress("%u sigs found\n", Count);

	if (OutputFileName == "")
		return;

	FILE *f = CreateStdioFile(OutputFileName);
	Progress("Writing %s\n", OutputFileName.c_str());
	WriteStdioFile(f, Buffer, FileSize);
	}
