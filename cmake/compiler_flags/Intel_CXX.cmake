# Special defines for Intel c++ compiler

set( ${PROJECT_NAME}_CXX_FLAGS_RELEASE "-O3 -xHost" )
set( ${PROJECT_NAME}_CXX_FLAGS_CLUSTER "-ipo" )

# Warning flags ...

set( ${PROJECT_NAME}_CXX_WARNING_FLAGS "-w3 -diag-disable:remark" )

# Treat 32bit/64bit compilation differently

set( ${PROJECT_NAME}_32BIT_CXX_FLAGS "-m32" )

set( ${PROJECT_NAME}_64BIT_CXX_FLAGS "-m64" )

