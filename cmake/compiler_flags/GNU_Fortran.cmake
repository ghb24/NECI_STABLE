# Special defines for gnu fortran compiler

set( ${PROJECT_NAME}_Fortran_FLAGS "-g -ffree-line-length-none -fPIC" )

set( ${PROJECT_NAME}_Fortran_FLAGS_DEBUG "-O0 -fbounds-check -fcheck=all -fbacktrace -finit-real=nan -ffpe-trap=invalid,zero,overflow,underflow" )
set( ${PROJECT_NAME}_Fortran_FLAGS_FASTDEBUG "-O2 -fbounds-check -fcheck=all -fbacktrace -finit-real=nan -ffpe-trap=invalid,zero,overflow,underflow" )
set( ${PROJECT_NAME}_Fortran_FLAGS_RELEASE "-O3 -march=native -mtune=native -funroll-loops" )
set( ${PROJECT_NAME}_Fortran_FLAGS_CLUSTER "-flto" )
set( ${PROJECT_NAME}_Fortran_LINKER_FLAGS_DEBUG "-rdynamic" )
set( ${PROJECT_NAME}_Fortran_LINKER_FLAGS_FASTDEBUG "-rdynamic" )

# Linker flags

set( ${PROJECT_NAME}_Fortran_LINKER_FLAGS_CLUSTER "-flto" )

# Warning flags. They are only set when warnings are enabled.
set( ${PROJECT_NAME}_Fortran_WARNING_FLAGS "-Wall -Wextra -Wno-zerotrip -Wno-maybe-uninitialized -Wno-do-subscript -Wcharacter-truncation -Wuse-without-only -Wimplicit-procedure -Winteger-division -Wreal-q-constant ")
# Disable some warnings for legacy code. **These flags are always enabled for the respective files.**.
set( ${PROJECT_NAME}_Fortran_relaxed_WARNING_FLAGS "-Wno-unused -Wno-character-truncation -Wno-use-without-only -Wno-implicit-procedure -Wno-integer-division -Wno-real-q-constant ")
if( CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10 )
    string(APPEND ${PROJECT_NAME}_Fortran_relaxed_WARNING_FLAGS " -fallow-argument-mismatch")
endif()

# Treat errors as warnings
set( ${PROJECT_NAME}_Fortran_WARN_ERROR_FLAG "-Werror")
# set( ${PROJECT_NAME}_Fortran_WARN_ERROR_FLAG "")

# Treat 32bit/64bit compilation differently

set( ${PROJECT_NAME}_32BIT_Fortran_FLAGS "-m32" )

set( ${PROJECT_NAME}_64BIT_Fortran_FLAGS "-m64" )

