# Special configuration to deal with oddities on Archer

set(CMAKE_SYSTEM_NAME Linux)

# Set the compiler. We don't need to use the force-detection as:
#
#  i) We are wanting to use the compiler wrappers (ftn), but we want the user to be able to supply to the
#     Archer build system which compilers that should be.
# ii) We don't want to add dependencies on specific versions of compilers and their linking needs

set( CMAKE_Fortran_COMPILER ftn )
set( CMAKE_CXX_COMPILER g++ )

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
# n.b. FORCE_ in the toolchain overrides the cmake/compiler_flags settings.
#      NECI_ does not overrider compiler_flags settings, but works with 

# Arbitrary compile flags can be overriden here.
# e.g.
#
# set( FORCE_CXX_FLAGS_RELEASE -O3 )   # Release mode specific CXX flags (also: DEBUG, CLUSTER)
# set( FORCE_Fortran_FLAGS -O3 )       # All build-type fortran flags

# Link flags

# Arbitrary compile flags can be overriden here.
# e.g.
#
# set( FORCE_CXX_LINK_FLAGS -O3 )

# Linker libraries

# Arbitrary flags are set here. Note that these are set by LANGUAGE and TARGET TYPE. They are set inside
# add_library or add_executable, and so only take impact if the LINKER_LANGUAGE setting is set on those
# types.
#
# e.g.
# set ( NECI_Fortran_STATIC_LINK_LIBRARIES fftw3 )     # Impacts static libraries
# set ( NECI_Fortran_EXE_LINK_LIBRARIES stdc++ )       # Impacts executables
# set ( NECI_Fortran_LINK_LIBRARIES m )                # Impacts both libraries and executables

# We don't want to propagate the g++ link elements through into the fortran linker.
# (these are found during the "Detecting CXX compiler ABI info" stage).
set( NECI_CXX_IMPLICIT_LINK_DIRECTORIES "" )



