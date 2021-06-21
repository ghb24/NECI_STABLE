# Special defines for gnu C++ compiler

set( ${PROJECT_NAME}_C_FLAGS "-fPIC" )
set( ${PROJECT_NAME}_C_FLAGS_CLUSTER "-flto" )

# Warning flags ...
set( ${PROJECT_NAME}_CXX_WARNING_FLAGS "-Wall -Wextra" )
# Treat errors as warnings
set( ${PROJECT_NAME}_CXX_WARN_ERROR_FLAG "-Werror")

# Treat 32bit/64bit compilation differently

set( ${PROJECT_NAME}_32BIT_CXX_FLAGS "-m32" )

set( ${PROJECT_NAME}_64BIT_CXX_FLAGS "-m64" )

