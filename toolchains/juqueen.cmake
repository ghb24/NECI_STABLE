set( CMAKE_SYSTEM_NAME Linux)
set( CMAKE_Fortran_COMPILER mpixlf90 )
set( CMAKE_CXX_COMPILER mpixlcxx )
set( CMAKE_C_COMPILER mpixlc )

set( NECI_FIND_LAPACK_NECI OFF )
set( LAPACK_NECI_LIBRARIES lapack esslbg )
set( LAPACK_NECI_FOUND true )
set( NECI_32BIT_Fortran_FLAGS "-q32 -qintsize=4 -qrealsize=4" )
set( NECI_64BIT_Fortran_FLAGS "-q64 -qintsize=8 -qrealsize=8" )
set( NECI_32BIT_C_FLAGS "-q32" )
set( NECI_64BIT_C_FLAGS "-q64" )
set( NECI_32BIT_CXX_FLAGS "-q32" )
set( NECI_64BIT_CXX_FLAGS "-q64" )
# Enforce 64bit architecture
set( CMAKE_SIZEOF_VOID_P 8 )

# Package overrides
set( NECI_FIND_MPI_NECI OFF )
set( MPI_NECI_LIBRARIES /bgsys/drivers/ppcfloor/comm/xl.ndebug )
#set( NECI_FIND_LAPACK_NECI OFF )
set( NECI_FIND_LibRT_NECI OFF )

# Compile flags
# =============
set( FORCE_CXX_FLAGS_RELEASE -O3 )
set( FORCE_Fortran_FLAGS -O3 )
set( FORCE_C_FLAGS -O3 )

#set( FORCE_CXX_FLAGS_RELEASE -O3 -qhot -qarch=qp -qtune=qp -lmass -lvmass -lmass_simd )
#set( FORCE_Fortran_FLAGS -O3 -qhot -qarch=qp -qtune=qp -qessl -lessl )
#set( FORCE_C_FLAGS -O3 -qhot -qarch=qp -qtune=qp -qessl -lessl )
# Link flags
# ==========
# Arbitrary compile flags can be overriden here.
set(  FORCE_CXX_LINK_FLAGS -O3 -Wl,--allow-multiple-definition)
set(  FORCE_Fortran_LINK_FLAGS -O3 -Wl,--allow-multiple-definition)

# Linker libraries
set( NECI_Fortran_STATIC_LINK_LIBRARIES stdc++ )

