# Special defines for Intel fortran compiler

set( ${PROJECT_NAME}_Fortran_FLAGS_DEBUG "-g -O0 -check bound" )
set( ${PROJECT_NAME}_Fortran_FLAGS_CLUSTER "-ipo" )

# Warning flags ...


# Treat 32bit/64bit compilation differently

set( ${PROJECT_NAME}_32BIT_Fortran_FLAGS "-m32" )

set( ${PROJECT_NAME}_64BIT_Fortran_FLAGS "-m64" )

