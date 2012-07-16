module dSFMT_interface

    ! A wrapper around a subset of the functionality in dSFMT RNG (double
    ! precision SIMD-oriented Fast Mersenne Twister Random Number Generator).
    !
    ! See also dSFMT.h/c


    use constants
    use iso_c_hack
    implicit none

    ! It is much faster to generate random numbers in blocks. genrand_real2 
    ! is a wrapper around accessing the random_store, filling it up again as 
    ! necessary.

    ! Testing indicates that 50000 is a very good size for the array storing 
    ! the random numbers.  Testing was done standalone, so undoubtedly !
    ! influenced by cache size and this might be different for real-world 
    ! applications, but it's easy to change to allocatable later on.

    integer, parameter :: random_store_size = 5 * 10**4
    real(c_double), save :: random_store(random_store_size)

    ! The next unused element in the store of random numbers.
    ! WARNING: random_store should be accessed via genrand_real2_dSFMT!
    integer, save :: current_element

!    real(dp), external :: genrand_close_open ! Given in dSFTM_wrapper.cpp.

    interface
        subroutine init_gen_rand_fwrapper (sd) bind(c)
            use iso_c_hack
            implicit none
            integer(c_int32_t), intent(in), value :: sd
        end subroutine
        subroutine fill_array_close_open_fwrapper(ar, sz) bind(c)
            use iso_c_hack
            implicit none
            real(c_double), intent(inout) :: ar(*)
            integer(c_int), intent(in), value :: sz
        end subroutine
    end interface


contains

    subroutine dSFMT_init(seed)

        ! Initialise the dSFMT RNG and fill random_store with
        ! a block of random numbers in interval [0,1).
        !
        ! In:
        !    seed: seed for the RNG.

        integer, intent(in) :: seed

        call init_gen_rand_fwrapper(int(seed,c_int32_t))

        call fill_array_close_open_fwrapper(random_store, &
                                            int(random_store_size, c_int))

        current_element = 1

    end subroutine dSFMT_init

    function genrand_real2_dSFMT() result(r)

        ! Return:
        !    random number in interval [0,1).  Name comes from the function
        !    defined in the original Mersenne Twist implementation.

        real(dp) :: r

        if (current_element == random_store_size+1) then
            ! Run out of random numbers: get more.
            current_element = 1
            call fill_array_close_open_fwrapper(random_store, &
                                                int(random_store_size, c_int))
        end if

        r = random_store(current_element)
        current_element = current_element + 1 

    end function genrand_real2_dSFMT

end module dSFMT_interface
