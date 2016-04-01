# Special defines for Intel C compiler

set( ${PROJECT_NAME}_C_FLAGS_CLUSTER "-ipo" )

# Warning flags ...

set( ${PROJECT_NAME}_C_WARNING_FLAGS "-warn all" )

# Treat 32bit/64bit compilation differently

set( ${PROJECT_NAME}_32BIT_C_FLAGS "-m32" )

set( ${PROJECT_NAME}_64BIT_C_FLAGS "-m64" )

