# Special defines for Intel c++ compiler

set( ${PROJECT_NAME}_CXX_FLAGS "-fPIC" )
set( ${PROJECT_NAME}_CXX_FLAGS_RELEASE "-O3 -xHost" )
set( ${PROJECT_NAME}_CXX_FLAGS_CLUSTER "-ipo" )

# Warning flags ...

# https://stackoverflow.com/questions/27343360/how-to-turn-on-icc-icpc-warnings
set( ${PROJECT_NAME}_CXX_WARNING_FLAGS "-Wall -Warray-bounds -Wchar-subscripts -Wcomment -Wenum-compare -Wformat -Wuninitialized -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull -Wparentheses -Wpointer-sign -Wreorder -Wreturn-type -Wsign-compare -Wsequence-point -Wtrigraphs -Wunused-function -Wunused-but-set-variable -Wunused-variable -Wwrite-strings -diag-disable:remark" )
# Treat errors as warnings
set( ${PROJECT_NAME}_CXX_WARN_ERROR_FLAG "-Werror")

# Treat 32bit/64bit compilation differently

set( ${PROJECT_NAME}_32BIT_CXX_FLAGS "-m32" )

set( ${PROJECT_NAME}_64BIT_CXX_FLAGS "-m64" )

