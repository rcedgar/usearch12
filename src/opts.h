#pragma once

#include "opt_enum.h"
#include <string>

#define PROGRAM_NAME	"usearch"
#define MY_VERSION	"12.0"

enum OTYPE
	{
	OTYPE_str,
	OTYPE_flt,
	OTYPE_uns,
	OTYPE_flag
	};

const string &oget_str(OENUM oe);
const char *oget_cstr(OENUM oe);

bool oget_flag(OENUM oe);

const string &oget_strd(OENUM oe, const string &dflt);
unsigned oget_unsd(OENUM oe, unsigned dflt);
double oget_fltd(OENUM oe, double dflt);

void oset_strd(OENUM oe, const string &dflt);
void oset_unsd(OENUM oe, unsigned dflt);
void oset_fltd(OENUM oe, double dflt);
void oset_flag(OENUM oe);

bool ofilled(OENUM oe);
bool ocmdline(OENUM oe);

unsigned oget_uns_(OENUM oe, const char *FN, uint LineNr);
double oget_flt_(OENUM oe, const char *FN, uint LineNr);
#define oget_uns(oe)	oget_uns_((oe), __FILE__, __LINE__)
#define oget_flt(oe)	oget_flt_((oe), __FILE__, __LINE__)