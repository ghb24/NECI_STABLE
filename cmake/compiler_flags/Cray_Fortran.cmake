# Special defines for Cray fortran compiler

set( ${PROJECT_NAME}_Fortran_FLAGS "-N255 -em" )
set( ${PROJECT_NAME}_Fortran_FLAGS_DEBUG "-R bcdps" )
set( ${PROJECT_NAME}_Fortran_FLAGS_FASTDEBUG "${${PROJECT_NAME}_Fortran_FLAGS_DEBUG}" )

# Treat 32bit/64bit compilation differently

set( ${PROJECT_NAME}_32BIT_Fortran_FLAGS "-s integer64 -s real64" )

set( ${PROJECT_NAME}_64BIT_Fortran_FLAGS "" )

