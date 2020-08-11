# Special configuration to deal with oddities on Maui

set(CMAKE_SYSTEM_NAME Linux)

set( FORCE_CXX_FLAGS_RELEASE "-O3 -xHost -std=c++11" )

set( CMAKE_Fortran_COMPILER ftn )
set( CMAKE_C_COMPILER cc )
set( CMAKE_CXX_COMPILER CC )

#Force the use of compiler wrappers also instead of mpiifort etc. This does the trick, apparently
set( MPI_Fortran_COMPILER ftn )
set( MPI_C_COMPILER cc )
set( MPI_CXX_COMPILER CC )


#This is supposed to prevent the FIND_MPI script from forcing the use of mpiifort etc.
set( ENABLE_INTEL_MPI_OVERRIDE OFF )

# Not sure if these actually do anything

set( NECI_FIND_MPI_NECI OFF )
set( NECI_FIND_LAPACK_NECI OFF )
set( NECI_FIND_LibRT_NECI OFF )

#set( NECI_SIMPLE_MPI_SEARCH OFF)

#set( MPI_NECI_LIBRARIES /opt/cray/pe/mpt/7.7.6/gni/mpich-intel/16.0/lib/libmpichf90_intel.so;/opt/cray/pe/mpt/7.7.6/gni/mpich-intel/16.0/lib/libfmpich_intel.so;/opt/cray/pe/mpt/7.7.6/gni/mpich-intel/16.0/lib/libmpichcxx_intel.so;/opt/cray/pe/mpt/7.7.6/gni/mpich-intel/16.0/lib/libmpl.so;/opt/cray/pe/mpt/7.7.6/gni/mpich-intel/16.0/lib/libtvmpich.so )

#set( CMAKE_FIND_LIBRARY_PREFIXES lib )
#set( CMAKE_FIND_LIBRARY_SUFFIXES .so .a )
#find_library( MPI_LIBRARIES mpich PATHS $MPICH_DIR/lib )
#find_library( MPI_LIBRARIES mpich PATHS /opt/cray/pe/mpt/7.7.6/gni/mpich-intel/16.0/lib )

#set( MPI_NECI_LIBRARIES "" )
#set( MPI_NECI_LIBRARIES ${MPI_LIBRARIES} )
#set( MPI_NECI_INCLUDE_PATH $MPICH_DIR/include )
#set( MPI_NECI_INCLUDE_PATH /opt/cray/pe/mpt/7.7.6/gni/mpich-intel/16.0/include )

#link_libraries( ${MPI_LIBRARIES} )
#target_link_libraries( ${MPI_LIBRARIES} )
#set( MPI_NECI_DEFINITIONS "" )
#set( MPI_NECI_FOUND true )


#set(HAVE_MPI ON)
#set( NECI_FIND_HDF5_NECI OFF )
#set( HDF5_NECI_LIBRARIES $HDF5_DIR/lib/libhdf5.a /$HDF5_DIR/lib/libhdf5_fortran.a )
#set( HDF5_NECI_DEFINITIONS "" )
#set( HDF5_NECI_INCLUDE_PATH $HDF5_DIR/include )
#set( HDF5_NECI_FOUND true )

