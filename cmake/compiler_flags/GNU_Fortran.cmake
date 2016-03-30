# Special defines for gnu fortran compiler

# Treat 32bit/64bit compilation differently

set( ${PROJECT_NAME}_32BIT_Fortran_FLAGS "-m32" )

set( ${PROJECT_NAME}_64BIT_Fortran_FLAGS "-m64 -fdefault-real-8 -fdefault-integer-8 -fdefault-double-8" )

