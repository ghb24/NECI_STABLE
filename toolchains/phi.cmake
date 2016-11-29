set(CMAKE_SYSTEM_NAME Linux)

set( CMAKE_Fortran_COMPILER mpiifort )
set( CMAKE_CXX_COMPILER mpiicpc )
set( CMAKE_C_COMPILER mpiicc )

set( NECI_32BIT_Fortran_FLAGS "-O3 -ipo -q32 -qintsize=4 -qrealsize=4 -xHost " )
set( NECI_64BIT_Fortran_FLAGS "-O3 -ipo -q64 -qintsize=8 -qrealsize=8 -xHost " )
# Enforce 64bit architecture
set( CMAKE_SIZEOF_VOID_P 8 )

set( FORCE_CXX_FLAGS_RELEASE "-O3 -xHost" )
set( FORCE_Fortran_FLAGS_RELEASE "-O3 -ipo -align array64byte -xHost " ) 
set( FORCE_C_FLAGS_RELEASE "-O3 -xHost" ) 
