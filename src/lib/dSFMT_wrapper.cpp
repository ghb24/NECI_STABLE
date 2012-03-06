#include "dSFMT.h"
#include "dSFMT_wrapper.h"

/** dsfmt internal state vector */
dsfmt_t dsfmt_global_data;

// Wrap around the required dSFMT functions so that they're accessible from 
// fortran. We use C++'s handy reference function to allow Fortran and C to 
// communicate, despite the different approaches in passing arguments.
// We only expose functions as needed.
//
// JSS (with a hat-tip to discussions with AJWT)

// This creates instantiations which may be linked, where before the
// declarations are inline static, and so don't appear in the object files

extern "C"
{
    void init_gen_rand_fwrapper(uint32_t seed)
    {
        // Initialise random number generator.
        // See also the main dSFMT code for the ability to initialise using an
        // array of seeds.
        init_gen_rand(seed);
    }

    double genrand_close_open_fwrapper(void)
    {
        // Return a random number in the interval [0,1).
        //
        // This uses the global state function.
        // dSFMT also has the ability to have "local" states---useful for
        // threading?
        return genrand_close_open();
    }

    void fill_array_close_open_fwrapper(double array[], int size)
    {
        // Fill an array of length size with random numbers in the interval
        // [0,1).
        // This is much faster than than repeated calls to genrand_close_open_.
        //
        // This uses the global state function.
        // dSFMT also has the ability to have "local" states---useful for
        // threading?
        fill_array_close_open(array, size);
    }

}
