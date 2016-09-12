# Special defines for the PGI C compiler

set( ${PROJECT_NAME}_C_FLAGS "-mcmodel=medium" )
set( ${PROJECT_NAME}_C_FLAGS_CLUSTER "-Mipa=fast" )

# Warning flags

set( ${PROJECT_NAME}_C_WARNING_FLAGS "-Minform=warn" )

# Treat 32bit/64bit compilation differently

set( ${PROJECT_NAME}_32BIT_C_FLAGS "-pc=32" )

set( ${PROJECT_NAME}_64BIT_C_FLAGS "-pc=64" )
