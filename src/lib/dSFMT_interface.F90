module dSFMT_interface

! This contains a handy wrapper around a subset of of the functionality in
! dSFMT RNG (double precision SIMD-oriented Fast Mersenne Twister random number
! generator).

! See also the functions defined in dSFMT_wrapper.cpp.

use constants
implicit none

! It is much faster to generate random numbers in blocks.
! genrand_real2 is a wrapper around accessing the random_store,
! filling it up again as necessary.

! Testing indicates that 50000 is a very good size for the array storing the
! random numbers.  Testing was done standalone, so undoubtedly influenced by
! cache size and this might be different for real-world applications, but it's easy to
! change to allocatable later on.
integer, parameter :: random_store_size=5*10**4
real(dp), save :: random_store(random_store_size)

integer, save :: current_element=1

real(dp), external :: genrand_close_open ! Given in dSFTM_wrapper.cpp.

contains

    subroutine dSFMT_init(seed)

        ! Initialise the dSFMT RNG and fill random_store with
        ! a block of random numbers in interval [0,1).
        !
        ! In:
        !    seed: seed for the RNG.

        integer, intent(in) :: seed
        integer :: ierr

        call init_gen_rand(seed)

        call fill_array_close_open(random_store, random_store_size)

    end subroutine dSFMT_init

    function genrand_real2_dSFMT() result(r)

        ! Return:
        !    random number in interval [0,1).  Name comes from the function
        !    defined in the original Mersenne Twist implementation.

        real(dp) :: r

        if (current_element == random_store_size+1) then
            ! Run out of random numbers: get more.
            current_element = 1
            call fill_array_close_open(random_store, random_store_size)
        end if

        r = random_store(current_element)
        current_element = current_element + 1 

    end function genrand_real2_dSFMT

    subroutine test_mt()

        use mt95
        real(4) :: t1(2), t2(2), s, etime
        real(dp) :: r
        integer :: i

        call genrand_init(7)
        s = etime(t1)
        do i = 1, 10**9
            call genrand_real2(r)
        end do
        s = etime(t2)
        write (6,*) 'mt95',r,t2-t1

        call init_gen_rand(7)
        s = etime(t1)
        do i = 1, 10**9
            r = genrand_close_open()
        end do
        s = etime(t2)
        write (6,*) 'dSFMT',r,t2-t1

        call dSFMT_init(7)
        s = etime(t1)
        do i = 1, 10**9
            r = genrand_real2_dSFMT()
        end do
        s = etime(t2)
        write (6,*) 'dSFMT2',r,t2-t1

        stop

    end subroutine test_mt

end module dSFMT_interface
