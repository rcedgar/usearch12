#ifndef onkmer_h
#define onkmer_h

typedef void fn_OnKmer(uint32 Code, uint32 Pos, void *UserData);
typedef void fn_OnKmerBoth(uint32 Code, uint32 Pos, bool Plus, void *UserData);

typedef uint fn_OnKmer_Skip(uint32 Code, uint32 Pos, void *UserData);
typedef uint fn_OnKmerBoth_Skip(uint32 Code, uint32 Pos, bool Plus, void *UserData);

#endif
