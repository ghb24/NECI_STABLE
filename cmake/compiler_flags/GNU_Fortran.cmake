# Special defines for gnu fortran compiler

set( ${PROJECT_NAME}_Fortran_FLAGS "-ffree-line-length-none -fPIC" )
# set( ${PROJECT_NAME}_Fortran_FLAGS "-g -ffree-line-length-none" )
set( ${PROJECT_NAME}_Fortran_FLAGS_DEBUG "-Og -fbounds-check -fcheck=all -fbacktrace -finit-real=nan -ffpe-trap=invalid,zero,overflow,underflow" )
set( ${PROJECT_NAME}_Fortran_FLAGS_RELEASE "-O3 -march=native -mtune=native -funroll-loops" )
set( ${PROJECT_NAME}_Fortran_FLAGS_CLUSTER "-flto" )
set( ${PROJECT_NAME}_Fortran_LINKER_FLAGS_DEBUG "-rdynamic" )

# Linker flags

set( ${PROJECT_NAME}_Fortran_LINKER_FLAGS_CLUSTER "-flto" )

# Warning flags ...
set( ${PROJECT_NAME}_Fortran_WARNING_FLAGS "-Wall -Wextra -Wno-zerotrip -Wno-maybe-uninitialized")
set( ${PROJECT_NAME}_Fortran_relaxed_WARNING_FLAGS "-Wno-unused -Wno-do-subscript")

# Treat errors as warnings
set( ${PROJECT_NAME}_Fortran_WARN_ERROR_FLAG "-Werror")

# Treat 32bit/64bit compilation differently

set( ${PROJECT_NAME}_32BIT_Fortran_FLAGS "-m32" )

set( ${PROJECT_NAME}_64BIT_Fortran_FLAGS "-m64" )

