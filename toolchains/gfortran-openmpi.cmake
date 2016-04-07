set(CMAKE_SYSTEM_NAME Linux)

# Enforce 64bit architecture

set( CMAKE_SIZEOF_VOID_P 8 )

# Compiler

include(CMakeForceCompiler)

CMAKE_FORCE_C_COMPILER       ( gcc GNU )
CMAKE_FORCE_CXX_COMPILER     ( g++ GNU )
CMAKE_FORCE_Fortran_COMPILER ( mpif90 GNU )


# Package overrides
# =================

# n.b. We can also define:
#
# set( ${PACKAGE_NAME}_FOUND ON )
# set( ${PACKAGE_NAME}_LIBRARIES ... )

set( NECI_FIND_MPI_NECI OFF )
set( NECI_FIND_LibRT OFF )

set( LIBRT_LIBRARIES rt )

# Compile flags
# =============

# n.b. FORCE_ in the toolchain overrides the cmake/compiler_flags settings.

# Arbitrary compile flags can be overriden here.
# e.g.
#
# set( FORCE_CXX_FLAGS_RELEASE -O3 )   # Release mode specific CXX flags (also: DEBUG, CLUSTER)
# set( FORCE_Fortran_FLAGS -O3 )       # All build-type fortran flags

# Link flags
# ==========

# Arbitrary compile flags can be overriden here.
# e.g.
#
# set(  FORCE_CXX_LINK_FLAGS -O3 )

# Linker libraries

# Arbitrary flags are set here. Note that these are set by LANGUAGE and TARGET TYPE. They are set inside
# add_library or add_executable, and so only take impact if the LINKER_LANGUAGE setting is set on those
# types.
#
# e.g.
# set ( NECI_Fortran_STATIC_LINK_LIBRARIES fftw3 )     # Impacts static libraries
# set ( NECI_Fortran_EXE_LINK_LIBRARIES stdc++ )       # Impacts executables
# set ( NECI_Fortran_LINK_LIBRARIES m )                # Impacts both libraries and executables

# By using CMakeForceCompiler, we break the autodetction of the c++ library requirement
# (STATIC --> adds these only to libneci, libkneci, ... not directly to the executables)
set( NECI_Fortran_STATIC_LINK_LIBRARIES stdc++ )
