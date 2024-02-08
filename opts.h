#pragma once

#include "opt_enum.h"
#include <string>

using namespace std;

//#define oget_str(x)		(string(""))
//#define oget_uns(x)		(0)
//#define oget_float(x)	(0.0)
//#define oget_flag(x)	(true)
//
//#define oget_strd(x, dflt)		(string(""))
//#define oget_unsd(x, dflt)		(0)
//#define oget_fltd(x, dflt)		(0.0)
//
//#define oset_flag(x)		/* empty */
//#define oset_int(x, dflt)	/* empty */
//#define oset_str(x, dflt)	/* empty */
//#define oset_flt(x, dflt)	/* empty */

const string &oget_str(OPT_ENUM oe);
const char *oget_cstr(OPT_ENUM oe);
unsigned oget_uns(OPT_ENUM oe);
double oget_flt(OPT_ENUM oe);
bool oget_flag(OPT_ENUM oe);

const string &oget_strd(OPT_ENUM oe, const string &dflt);
unsigned oget_unsd(OPT_ENUM oe, unsigned dflt);
double oget_fltd(OPT_ENUM oe, double dflt);

void oset_uns(OPT_ENUM oe, unsigned value);

void oset_strd(OPT_ENUM oe, const string &dflt);
void oset_unsd(OPT_ENUM oe, unsigned dflt);
void oset_fltd(OPT_ENUM oe, double dflt);
void oset_flag(OPT_ENUM oe);

bool ofilled_str(OPT_ENUM oe);
bool ofilled_uns(OPT_ENUM oe);
bool ofilled_flt(OPT_ENUM oe);
bool ofilled_flag(OPT_ENUM oe);
