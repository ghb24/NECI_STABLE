# Special defines for Intel fortran compiler

set( ${PROJECT_NAME}_Fortran_FLAGS_DEBUG "-g -O0 -check bound" )
set( ${PROJECT_NAME}_Fortran_FLAGS_CLUSTER "-ipo" )

# Warning flags ...

set( ${PROJECT_NAME}_Fortran_WARNING_FLAGS "-warn all" )

# Linker flags

set( ${PROJECT_NAME}_Fortran_LINKER_FLAGS_CLUSTER "-ipo" )

# A particular ifort fix (I can't find any other reference to this CMake variable...)

# -i_dynamic is incorrectly added in CMake ifort configuration
set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")

# Treat 32bit/64bit compilation differently

set( ${PROJECT_NAME}_32BIT_Fortran_FLAGS "-m32" )

set( ${PROJECT_NAME}_64BIT_Fortran_FLAGS "-m64" )
