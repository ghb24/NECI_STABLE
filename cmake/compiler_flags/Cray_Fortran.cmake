# Special defines for Cray fortran compiler

# Treat 32bit/64bit compilation differently

set( ${PROJECT_NAME}_32BIT_Fortran_FLAGS "-s integer64 -s real64" )

set( ${PROJECT_NAME}_64BIT_Fortran_FLAGS "" )

