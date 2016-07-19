set( CMAKE_SYSTEM_NAME Linux)
set( CMAKE_Fortran_COMPILER mpixlf90_r )
set( CMAKE_CXX_COMPILER mpixlcxx_r )
set( CMAKE_C_COMPILER mpixlc_r )

set( NECI_FIND_LAPACK_NECI OFF )
set( LAPACK_NECI_FOUND true )
set( CMAKE_FIND_LIBRARY_PREFIXES lib )
set( CMAKE_FIND_LIBRARY_SUFFIXES .so .a )
find_library( ESSLBG_LIBRARY esslbg HINTS /bgsys/local/lib )
find_library( LAPACK_LIBRARY lapack HINTS /bgsys/local/lapack/3.6.0/lib )
set( LAPACK_NECI_LIBRARIES "" )
link_libraries( ${LAPACK_LIBRARY} ${ESSLBG_LIBRARY} )
#set( LAPACK_NECI_LIBRARIES ${ESSLBG_LIBRARY} ${LAPACK_LIBRARY} )
#target_link_libraries( ${LAPACK_LIBRARY} ${ESSLBG_LIBRARY} )
#This should fix the flush_ and the hostnm_ for IBM libraries
add_definitions( -DBLUEGENE_HACKS )

set( NECI_32BIT_Fortran_FLAGS "-q32 -qintsize=4 -qrealsize=4" )
set( NECI_64BIT_Fortran_FLAGS "-q64 -qintsize=8 -qrealsize=8" )
set( NECI_F77_FLAGS -qfixed )
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

set( FORCE_CXX_FLAGS_DEBUG "-O0 -g" )
set( FORCE_Fortran_FLAGS_DEBUG "-O0 -g -qdpcl -qcheck -qextchk -qdbg" )
set( FORCE_C_FLAGS_DEBUG "-O0 -g" )

set( FORCE_CXX_FLAGS_RELEASE -O3 -qhot -qarch=qp -qtune=qp -lmass -lvmass -lmass_simd )
set( FORCE_Fortran_FLAGS_RELEASE -O3 -qhot -qarch=qp -qtune=qp -qessl -lessl )
set( FORCE_C_FLAGS_RELEASE -O3 -qhot -qarch=qp -qtune=qp -qessl -lessl )

# Link flags

set(  FORCE_CXX_LINKER_FLAGS -Wl,--allow-multiple-definition)
set(  FORCE_Fortran_LINKER_FLAGS -Wl,--allow-multiple-definition)

# Linker libraries

set( NECI_Fortran_STATIC_LINK_LIBRARIES stdc++ )

set( NECI_DISABLE_SSE2 ON )
