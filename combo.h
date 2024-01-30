#ifndef combo_h
#define combo_h

typedef void (*ON_SET)(const vector<unsigned> &v);
typedef vector<unsigned>::const_iterator vit;

void EnumPowerSetPerms(unsigned N, ON_SET OnPerm);

#endif // combo_h
