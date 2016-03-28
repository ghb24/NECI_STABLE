set(CMAKE_SYSTEM_NAME Linux)

# Enforce 64bit architecture

set( CMAKE_SIZEOF_VOID_P 8 )

# Compiler

include(CMakeForceCompiler)

CMAKE_FORCE_C_COMPILER       ( gcc GNU )
CMAKE_FORCE_CXX_COMPILER     ( g++ GNU )
CMAKE_FORCE_Fortran_COMPILER ( mpif90 GNU )

# Set the default implicit linker flags for the different languages 

#set( CMAKE_IMPLICIT_CXX_LINK_LIBRARIES stdc++ m c )
#set( CMAKE_IMPLICIT_Fortran_LINK_LIBRARIES gfortran m quadmath m c )

# TODO: Compiler flags!!!

# Override MPI searching. Override arbitrary package searching...

set( NECI_FIND_MPI OFF )
set( NECI_FIND_LibRT OFF )

# Compile flags

set( NECI_CXX_FLAGS_RELEASE -O9 -help )

# Cluster compilation flags

# Debug flags

# Link flags

# By using CMakeForceCompiler, we break the autodetction of the c++ library requirement
set( NECI_Fortran_STATIC_LINK_LIBRARIES stdc++ rt )
