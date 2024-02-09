#pragma once

#include "opt_enum.h"
#include <string>

#define PROGRAM_NAME	"usearch"
#define MY_VERSION	"12.0"

enum OPT_TYPE
	{
	OPT_TYPE_str,
	OPT_TYPE_flt,
	OPT_TYPE_uns,
	OPT_TYPE_flag
	};

const string &oget_str(OPT_ENUM oe);
const char *oget_cstr(OPT_ENUM oe);
unsigned oget_uns(OPT_ENUM oe);
double oget_flt(OPT_ENUM oe);
bool oget_flag(OPT_ENUM oe);

const string &oget_strd(OPT_ENUM oe, const string &dflt);
unsigned oget_unsd(OPT_ENUM oe, unsigned dflt);
double oget_fltd(OPT_ENUM oe, double dflt);

void oset_strd(OPT_ENUM oe, const string &dflt);
void oset_unsd(OPT_ENUM oe, unsigned dflt);
void oset_fltd(OPT_ENUM oe, double dflt);
void oset_flag(OPT_ENUM oe);

bool ofilled(OPT_ENUM oe);
