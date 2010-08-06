//-----------------------------------------------------------------------------
// MurmurHash2, by Austin Appleby

// Note - This code makes a few assumptions about how your machine behaves -

// 1. We can read a 4-byte value from any address without crashing
// 2. sizeof(int) == 4

// And it has a few limitations -

// 1. It will not work incrementally.
// 2. It will not produce the same results on little-endian and big-endian
//    machines.

#include <stdio.h>

unsigned int MurmurHash2 ( const void * key, int len, unsigned int seed )
{
    // 'm' and 'r' are mixing constants generated offline.
    // They're not really 'magic', they just happen to work well.

    const unsigned int m = 0x5bd1e995;
    const int r = 24;

    // Initialize the hash to a 'random' value

    unsigned int h = seed ^ len;

    // Mix 4 bytes at a time into the hash

    const unsigned char * data = (const unsigned char *)key;

    while(len >= 4)
    {
        unsigned int k = *(unsigned int *)data;

        k *= m; 
        k ^= k >> r; 
        k *= m; 
        
        h *= m; 
        h ^= k;

        data += 4;
        len -= 4;
    }
    
    // Handle the last few bytes of the input array


    switch(len)
    {
    case 3: h ^= data[2] << 16;
    case 2: h ^= data[1] << 8;
    case 1: h ^= data[0];
            h *= m;
    };

    // Do a few final mixes of the hash to ensure the last few
    // bytes are well-incorporated.

    h ^= h >> 13;
    h *= m;
    h ^= h >> 15;

    return h;
} 
// Wrapper to be used from fortran.
// 1. fortran appends an underscore to the object names (C doesn't) so we have to add it explicitly.
// 2. fortran passes everything by reference by default.  C passes only arrays by reference by default.
//    MurmurHash takes the length of the array (key) and the seed as integers.
//    The wrapper takes the arguments as references and passes the references or their values
//    (as appropriate) to MurmurHash and returns the hash.
unsigned int murmurhash2wrapper_ ( int *key, int *len, unsigned int *seed ) {
    // In:
    //    key: array to be hashed.
    //    len: length of array (in units of 4 bytes)
    //    seed: random seed used to generate the hash.
    // Need to dereference len and seed.
    return MurmurHash2(key,*len,*seed);
}
