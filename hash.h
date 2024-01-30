#ifndef hash_h
#define hash_h

/***
Hash functions from:
https://stackoverflow.com/questions/664014/what-integer-hash-function-are-good-that-accepts-an-integer-hash-key

References:
http://xorshift.di.unimi.it/splitmix64.c
http://zimbry.blogspot.it/2011/09/better-bit-mixing-improving-on.html
***/

static inline uint64 Hash64(uint64 x)
    {
    x = (x ^ (x >> 30)) * uint64(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * uint64(0x94d049bb133111eb);
    x = x ^ (x >> 31);
    return x;
    }

/***
According to stackoverflow answer:
"The following algorithm provides a very good statistical distribution.
Each input bit affects each output bit with about 50% probability.
There are no collisions (each input results in a different output).

You can reverse the process (get the input value from the hash) if you replace the 0x45d9f3b with 0x119de1f3 (the multiplicative inverse):

unsigned int unhash(unsigned int x) {
    x = ((x >> 16) ^ x) * 0x119de1f3;
    x = ((x >> 16) ^ x) * 0x119de1f3;
    x = (x >> 16) ^ x;
    return x;
}
"
***/

static inline uint32 Hash32(uint32 x)
    {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
    }

#endif // hash_h
