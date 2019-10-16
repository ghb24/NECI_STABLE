# Special defines for gnu C compiler

set( ${PROJECT_NAME}_C_FLAGS_CLUSTER "-flto" )

# Warning flags ...
set( ${PROJECT_NAME}_C_WARNING_FLAGS "-Wall -Wextra" )
# Treat errors as warnings
set( ${PROJECT_NAME}_C_WARN_ERROR_FLAG "-Werror")

# Treat 32bit/64bit compilation differently

set( ${PROJECT_NAME}_32BIT_C_FLAGS "-m32" )

set( ${PROJECT_NAME}_64BIT_C_FLAGS "-m64" )

