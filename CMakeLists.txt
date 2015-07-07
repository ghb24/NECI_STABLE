cmake_minimum_required(VERSION 2.8.4)
project(neci Fortran CXX C)

# Add our custom Modules directory, so that we can search for libraries
# that are not provided with finders by default
# (see www.cmake.org/Wiki/CMake:How_To_Find_Libraries)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

# Consider the possible source files
file(GLOB_RECURSE F90_FILES src/*.F90 src/*.F src/*.cpp)
file(GLOB_RECURSE F90_TEMPLATE_FILES src/*.F90.template)
list(REMOVE_ITEM F90_FILES ${CMAKE_SOURCE_DIR}/src/libstub.F)

file(GLOB_RECURSE UNIT_TEST_FILES unit_tests/*.F90 unit_tests/*.f90)
set(UNIT_TEST_FILES ${UNIT_TEST_FILES} ${F90_FILES})
list(REMOVE_ITEM UNIT_TEST_FILES ${CMAKE_SOURCE_DIR}/src/necimain.F90)



##############################################################################
#                                                                            #
# Compiler specific configuration                                            #
#                                                                            #
##############################################################################

message("Configuring Fortran compilers")
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
	message("-- Using ifort")
	# -i_dynamic is incorrectly added in CMake ifort configuration
	set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")
	add_definitions(-D__IFORT)

	set (CMAKE_Fortran_FLAGS_DEBUG "-check bounds")

elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

	message("-- Using gfortran")
	add_definitions(-D__GFORTRAN__)

	set(CMAKE_Fortran_FLAGS_DEBUG "-fbounds-check")

else()

	message (FATAL_ERROR "Fortran compiler unknown: ${CMAKE_Fortran_COMPILER_ID}")

endif()

#
# Deal with the additional build step introduced by the templated
# files
set(TEMPLATE_TARGETS "")
foreach(template_file ${F90_TEMPLATE_FILES})

	# For each of the templated files, we create a custom build command
	# command which depends on the templated file, and produces a .F90
	# file which is added to the list of dependencies for the executable.
	get_filename_component(base_filename ${template_file} NAME_WE)
	set(target_filename ${CMAKE_BINARY_DIR}/templated/${base_filename}.F90)
	set(TEMPLATE_TARGETS ${TEMPLATE_TARGETS} ${target_filename})

	add_custom_command(
		OUTPUT ${target_filename}
		COMMAND ${CMAKE_SOURCE_DIR}/tools/f90_template.py ${template_file} ${target_filename}
		DEPENDS ${template_file}
	)

endforeach(template_file)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/templated)

#######################################################################
#                                                                     #
# Buildable targets                                                   #
#                                                                     #
#######################################################################

#
# Build the executable!
add_executable(neci ${F90_FILES} ${TEMPLATE_TARGETS})
add_executable(mneci ${F90_FILES} ${TEMPLATE_TARGETS})
add_executable(dneci ${F90_FILES} ${TEMPLATE_TARGETS})
add_executable(kneci ${F90_FILES} ${TEMPLATE_TARGETS})

# Add the testing executables
add_executable(neci_test ${TEMPLATE_TARGETS} ${UNIT_TEST_FILES})
add_executable(mneci_test ${TEMPLATE_TARGETS} ${UNIT_TEST_FILES})
add_executable(dneci_test ${TEMPLATE_TARGETS} ${UNIT_TEST_FILES})
add_executable(kneci_test ${TEMPLATE_TARGETS} ${UNIT_TEST_FILES})

# List of targets. Then loop over repeated target properties
set(TARGETS neci mneci dneci kneci neci_test mneci_test dneci_test kneci_test)
foreach(target ${TARGETS})

	# Set the fortran modules directory for cleanliness, and to avoid collisions
	# between targets
	set_target_properties(${target} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules/${target})

	# Ensure that we use the Fortran compiler for linking
	set_target_properties(${target} PROPERTIES LINKER_LANGUAGE Fortran)

endforeach(target)

#
# And add the target-specific compler flags.
target_compile_definitions(neci PUBLIC -DHElement_t=real\(dp\))
target_compile_definitions(dneci PUBLIC -DHElement_t=real\(dp\) -D__DOUBLERUN)
target_compile_definitions(mneci PUBLIC -DHElement_t=real\(dp\) -D__PROG_NUMRUNS)
target_compile_definitions(kneci PUBLIC -DHElement_t=complex\(dp\) -D__CMPLX)
target_compile_definitions(neci_test PUBLIC -DHElement_t=real\(dp\))
target_compile_definitions(dneci_test PUBLIC -DHElement_t=real\(dp\) -D__DOUBLERUN)
target_compile_definitions(mneci_test PUBLIC -DHElement_t=real\(dp\) -D__PROG_NUMRUNS)
target_compile_definitions(kneci_test PUBLIC -DHElement_t=complex\(dp\) -D__CMPLX)

# TODO:
# Test for compiler features, to enable specific flags
# Add librt, libm as required
# Enable shared memory
# Select optimised, debug or production builds


#
# If we can get the repository details, then do so.
message("VCS configuration")
find_package(Git)
if (${Git_FOUND})
	message("-- Git has been found")
	Git_WC_INFO("." "DIR")
	add_definitions(-D_VCS_VER='${DIR_WC_REVISION_HASH}')
	message("-- SHAID: ${DIR_WC_REVISION_HASH}")
	if (${DIR_WC_CHANGES} OR ${DIR_WC_CACHED_CHANGES})
		add_definitions(-D_WORKING_DIR_CHANGES)
		message("-- Working directory contains changes")
	endif()
else()
	add_definitions(-D_VCS_VER='[Unknown]')
	message("-- Could NOT find git")
endif()

# Add the compile definitions
add_definitions(-D__Linux)
add_definitions(-DPOINTER8)
add_definitions(-DHAVE_SSE2)

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
	add_definitions(-D__INT64)
	add_definitions(-DPOINTER8)
endif()


add_definitions(-DMAXMEM=99999)
add_definitions(-DDSFMT_MEXP=19937)
add_definitions(-D_CONFIG="CMake-new")


include_directories(BEFORE SYSTEM ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/unit_tests)

#
# Some form of lapack is always needed (be that acml, mkl, LAPACK, ...)
message("Configuring BLAS/LAPACK")
find_package(LAPACK REQUIRED)
set(NECI_LIBRARIES ${LAPACK_LIBRARIES})
message("-- Using libraries: ${LAPACK_LIBRARIES}")

# If FFTW is present, then we want to use it. Otherwise disable it in NECI
message("Configuring FFTW")
find_package(FFTW)
if (${FFTW_FOUND})
	message("-- Using libraries: ${FFTW_LIBRARIES}")
	set(NECI_LIBRARIES ${NECI_LIBRARIES} ${FFTW_LIBRARIES})
else()
	message("-- Disabling FFTW")
	add_definitions(-DDISABLE_FFTW)
endif()

# Look for MPI...
find_package(MPI)
set(NECI_LIBRARIES ${NECI_LIBRARIES} ${MPI_Fortran_LIBRARIES})

# Ensure that all of our executables link correctly
foreach(target ${TARGETS})
	target_link_libraries(${target} ${NECI_LIBRARIES})
endforeach()
