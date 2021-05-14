
# Get the python interpreter
set(f2py_compiler ${MPI_Fortran_COMPILER})
get_target_property(NECI_MODULE_DIR libneci Fortran_MODULE_DIRECTORY)
get_target_property(NECI_LIB_DIR libneci LIBRARY_OUTPUT_DIRECTORY)
# Setup the script
configure_file(invoke_f2py.sh.in invoke_f2py.sh @ONLY)

add_custom_target(neci_guga_pylib
  COMMAND bash invoke_f2py.sh
  DEPENDS libneci
  BYPRODUCTS neci_guga.so
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  VERBATIM
)
