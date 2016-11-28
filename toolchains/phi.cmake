set(CMAKE_SYSTEM_NAME Linux)

set( CMAKE_Fortran_COMPILER mpiifort )
set( CMAKE_CXX_COMPILER mpiicpc )
set( CMAKE_C_COMPILER mpiicc )

set( NECI_32BIT_Fortran_FLAGS "-O3 -xHost -q32 -qintsize=4 -qrealsize=4" )
set( NECI_64BIT_Fortran_FLAGS "-O3 -xHost -q64 -qintsize=8 -qrealsize=8" )
# Enforce 64bit architecture
set( CMAKE_SIZEOF_VOID_P 8 )

set( FORCE_CXX_FLAGS_RELEASE "-O3 -xHost" )
set( FORCE_Fortran_FLAGS_RELEASE "-O3 -xHost -align array64byte" ) 
set( FORCE_C_FLAGS_RELEASE "-O3 -xHost" ) 
