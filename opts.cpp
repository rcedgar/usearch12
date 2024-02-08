#include "myutils.h"

const char *o_to_str(OPT_ENUM o)
	{
	switch (o)
		{
#define o(x)	case OPT_##x: return #x;
#include "o_all.h"
		}
	return "OPT_?";
	}

const string &oget_str(OPT_ENUM oe)
	{
	static string str;
	return str;
	}

unsigned oget_uns(OPT_ENUM oe)
	{
	return 0;
	}

double oget_flt(OPT_ENUM oe)
	{
	return 0;
	}

bool oget_flag(OPT_ENUM oe)
	{
	return false;
	}

const string &oget_strd(OPT_ENUM oe, const string &dflt)
	{
	return dflt;
	}

unsigned oget_unsd(OPT_ENUM oe, unsigned dflt)
	{
	return dflt;
	}

double oget_fltd(OPT_ENUM oe, double dflt)
	{
	return dflt;
	}

void oset_strd(OPT_ENUM oe, const string &dflt)
	{
	}

void oset_unsd(OPT_ENUM oe, unsigned dflt)
	{
	}

void oset_uns(OPT_ENUM oe, unsigned dflt)
	{
	}

void oset_fltd(OPT_ENUM oe, double dflt)
	{
	}

void oset_flag(OPT_ENUM oe)
	{
	}

bool ofilled_str(OPT_ENUM oe)
	{
	return false;
	}

bool ofilled_uns(OPT_ENUM oe)
	{
	return false;
	}

bool ofilled_flt(OPT_ENUM oe)
	{
	return false;
	}

bool ofilled_flag(OPT_ENUM oe)
	{
	return false;
	}

const char *oget_cstr(OPT_ENUM oe)
	{
	return 0;
	}
