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

get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
message("CMake compiler " ${CMAKE_Fortran_COMPILER})
message("CMake binary_dir " ${CMAKE_BINARY_DIR})
message("CMake project_binary_dir " ${CMAKE_PROJECT_BINARY_DIR})
message("CMake build dir: " ${CMAKE_BUILD_DIRECTORY})
message("CMake build files dir: " ${CMAKE_BUILD_FILES_DIRECTORY})
message("CMake runtime output: " ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
message("CMake executable path: " ${EXECUTABLE_OUTPUT_PATH})
message("CMake cachefile dir: " ${CMAKE_CACHEFILE_DIR})
message("CMake compiler name: " ${CMAKE_Fortran_COMPILER})
message("CMake compiler name: " ${Fortran_COMPILER_NAME})
message("CMake compiler name: " ${CMAKE_C_COMPILER})
message("CMake compiler name: " ${CMAKE_CXX_COMPILER})
message("Shared link: " ${CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS})

# Deal with compiler specific issues
if (Fortran_COMPILER_NAME MATCHES "ifort.*")
	# -i_dynamic is incorrectly added in CMake ifort configuration
	set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")
endif()
message("Shared link: " ${CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS})

#
# Deal with the additional build step introduced by the templated
# files
set(TEMPLATE_TARGETS "")
foreach(template_file ${F90_TEMPLATE_FILES})

	# For each of the templated files, we create a custom build command
	# command which depends on the templated file, and produces a .F90
	# file which is added to the list of dependencies for the executable.
	get_filename_component(base_filename ${template_file} NAME_WE)
	set(target_filename ${base_filename}.F90)
	set(TEMPLATE_TARGETS ${TEMPLATE_TARGETS} ${target_filename})

	add_custom_command(
		OUTPUT ${target_filename}
		COMMAND ${CMAKE_SOURCE_DIR}/tools/f90_template.py ${template_file} ${target_filename}
		DEPENDS ${template_file}
	)

endforeach(template_file)

#
# Build the executable!
add_executable(neci ${F90_FILES} ${TEMPLATE_TARGETS})
add_executable(mneci ${F90_FILES} ${TEMPLATE_TARGETS})
add_executable(dneci ${F90_FILES} ${TEMPLATE_TARGETS})
add_executable(kneci ${F90_FILES} ${TEMPLATE_TARGETS})

# Neci is fundamentally a fortran program, so don't default to using
# C to link things...
set_target_properties(neci PROPERTIES LINKER_LANGUAGE Fortran)
set_target_properties(mneci PROPERTIES LINKER_LANGUAGE Fortran)
set_target_properties(dneci PROPERTIES LINKER_LANGUAGE Fortran)
set_target_properties(kneci PROPERTIES LINKER_LANGUAGE Fortran)

# TODO:
# Test for compiler features, to enable specific flags
# Add librt, libm as required
# Enable shared memory
# Be able to build kneci, dneci., mneci
# Select optimised, debug or production builds


#
# If we can get the repository details, then do so.
find_package(Git)
if (${Git_FOUND})
	Git_WC_INFO("." "DIR")
	add_definitions(-D_VCS_VER='${DIR_WC_REVISION_HASH}')
	message("-- Git has been found")
	message("-- Git revision: ${DIR_WC_REVISION_HASH}")
	message("-- Git name: ${DIR_WC_REVISION_NAME}")
	message("-- Git changes: ${DIR_WC_CHANGES} ${DIR_WC_CACHED_CHANGES}")
else()
	add_definitions(-D_VCS_VER='[Unknown]')
	message("-- Could NOT find git")
endif()

# Add the compile definitions
add_definitions(-D__Linux)
add_definitions(-DPOINTER8)
add_definitions(-DHAVE_SSE2)
add_definitions(-D__INT64)
#add_definitions(-D__IFORT)
add_definitions(-D__GFORTRAN__)
add_definitions(-DMAXMEM=99999)
add_definitions(-DDSFMT_MEXP=19937)
add_definitions(-D_CONFIG="CMake-new")

# This is a real build. Consider kneci...
target_compile_definitions(neci PUBLIC -DHElement_t=real\(dp\))
target_compile_definitions(dneci PUBLIC -DHElement_t=real\(dp\))
target_compile_definitions(mneci PUBLIC -DHElement_t=real\(dp\))
target_compile_definitions(kneci PUBLIC -DHElement_t=complex\(dp\))

include_directories(BEFORE SYSTEM ${CMAKE_SOURCE_DIR}/src)

#
# Some form of lapack is always needed (be that acml, mkl, LAPACK, ...)
find_package(LAPACK REQUIRED)
set(NECI_LIBRARIES ${LAPACK_LIBRARIES})
message("Using LAPACK libraries: " ${LAPACK_LIBRARIES})

# If FFTW is present, then we want to use it. Otherwise disable it in NECI
find_package(FFTW)
if (${FFTW_FOUND})
	set(NECI_LIBRARIES ${NECI_LIBRARIES} ${FFTW_LIBRARIES})
else()
	message("Disabling FFTW")
	add_definitions(-DDISABLE_FFTW)
endif()

# Look for MPI...
find_package(MPI)
set(NECI_LIBRARIES ${NECI_LIBRARIES} ${MPI_Fortran_LIBRARIES})

# Ensure that all of our executables link correctly
target_link_libraries(neci ${NECI_LIBRARIES})
target_link_libraries(mneci ${NECI_LIBRARIES})
target_link_libraries(dneci ${NECI_LIBRARIES})
target_link_libraries(kneci ${NECI_LIBRARIES})
