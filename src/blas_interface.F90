#include "macros.h"

module blas_interface_mod
    !! This module should contain the declaration of BLAS routines.
    !!
    !! At the moment they are just declared as external,
    !! but one can change that one by one.
    better_implicit_none

    private

    public :: dcopy, dgemm, dsyev, zheev, zgemm, dgetrf, dgetri, dgeev, dscal

    external :: dcopy, dgemm, dsyev, zheev, zgemm, dgetrf, dgetri, dgeev, dscal

end module
