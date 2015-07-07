#ifndef dSFMT_wrapper_H
#define dSFMT_wrapper_H

// See notes in dSFMT_wrapper.cpp.

#ifdef __cplusplus
extern "C"
{
#endif

    void init_gen_rand_(uint32_t &seed);

    double genrand_close_open_(void);

    void fill_array_close_open_(double array[], int &size);

#ifdef __cplusplus
}
#endif

#endif // dSFMT_wrapper_H
