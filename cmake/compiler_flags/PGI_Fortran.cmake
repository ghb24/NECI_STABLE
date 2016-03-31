# Special defines for the PGI fortran compiler

set( ${PROJECT_NAME}_Fortran_FLAGS "-mcmodel=medium -Msignextend" )
set( ${PROJECT_NAME}_Fortran_FLAGS_DEBUG "-Mbounds" )
set( ${PROJECT_NAME}_Fortran_FLAGS_RELEASE "-tp x64 -fastsse" )

set( ${PROJECT_NAME}_Fortran_LINK_FLAGS "-mcmodel=medium" )

# Warning flags

set( ${PROJECT_NAME}_Fortran_WARNING_FLAGS "-Minform=warn" )

# Treat 32bit/64bit compilation differently

set( ${PROJECT_NAME}_32BIT_Fortran_FLAGS "-pc=32" )

set( ${PROJECT_NAME}_64BIT_Fortran_FLAGS "-pc=64 -r8 -i8" )

# The linker command -lpthread doesn't work for the PGI linker (due to some
# internal trickery with the pthread library
# --> Force the rest of the code to think that libpthreads has already been
#     found, and then add the magic -pthread command
# SDS: This no longer seems needed for pgfortran, which does the work automagically
# set(ALL_LINKER_FLAGS "${ALL_LINKER_FLAGS} -pthread")
