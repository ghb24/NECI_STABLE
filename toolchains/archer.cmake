# Special configuration to deal with oddities on Archer

set(CMAKE_SYSTEM_NAME Linux)

# We don't need to forcibly set the compiler. That will happen just fine automatically.


# Some packages are automagically included on archer.

set( NECI_FIND_MPI_NECI OFF )
set( NECI_FIND_LAPACK_NECI OFF )
set( NECI_FIND_LibRT_NECI OFF )


# Ensure that the loaded module is used for hdf5.
# Normally we would link against automatically discovered .so files, but hydra seems to mess this up
# as the module has been compiled using Position Independent Code.
set( NECI_FIND_HDF5_NECI OFF )
set( HDF5_NECI_LIBRARIES $ENV{HDF5_ROOT}/lib/libhdf5.a $ENV{HDF5_ROOT}/lib/libhdf5_fortran.a )
set( HDF5_NECI_DEFINITIONS "" )
set( HDF5_NECI_INCLUDE_PATH $ENV{HDF5_ROOT}/include )
set( HDF5_NECI_FOUND true )

# Compile flags

# Arbitrary compile flags can be overriden here.
# e.g.
#
# set( NECI_CXX_FLAGS_RELEASE -O3 )   # Release mode specific CXX flags (also: DEBUG, CLUSTER)
# set( NECI_Fortran_FLAGS -O3 )       # All build-type fortran flags

# Link flags

# Arbitrary compile flags can be overriden here.
# e.g.
#
# set( NECI_CXX_LINK_FLAGS -O3 )

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
set( NECI_Fortran_STATIC_LINK_LIBRARIES stdc++ rt )

