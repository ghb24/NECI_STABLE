#include "macros.h"

module blas_interface_mod
    better_implicit_none

    private

    public :: dcopy, dgemm, dsyev, zheev, zgemm, dgetrf, dgetri, dgeev

    external :: dcopy, dgemm, dsyev, zheev, zgemm, dgetrf, dgetri, dgeev

end module
