# Special defines for gnu fortran compiler

set( ${PROJECT_NAME}_Fortran_FLAGS "-ffree-line-length-none" )
set( ${PROJECT_NAME}_Fortran_FLAGS_DEBUG "-g -O0 -fbounds-check" )
set( ${PROJECT_NAME}_Fortran_FLAGS_RELEASE "-O3" )
set( ${PROJECT_NAME}_Fortran_FLAGS_CLUSTER "-flto" )

# Warning flags ...

set( ${PROJECT_NAME}_Fortran_WARNING_FLAGS "-Wall -Wextra" )

# Treat 32bit/64bit compilation differently

set( ${PROJECT_NAME}_32BIT_Fortran_FLAGS "-m32" )

set( ${PROJECT_NAME}_64BIT_Fortran_FLAGS "-m64 -fdefault-real-8 -fdefault-integer-8 -fdefault-double-8" )

