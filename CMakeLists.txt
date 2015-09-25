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
# User controllable options                                                  #
#                                                                            #
##############################################################################

# The build type, specifies the optimisation level
if (NOT DEFINED CMAKE_BUILD_TYPE OR NOT CMAKE_BUILD_TYPE)
	message("Setting type ${CMAKE_BUILD_TYPE}")
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING
		"Type of build, options are: Debug Release RelWithDebInfo Cluster" FORCE)
else()
	message("Defined type --${CMAKE_BUILD_TYPE}--")
endif()

option(PARALLEL "Enable parallel build (MPI)" on)
option(HDF "Enable HDF5 popsfiles" on)
option(WARNINGS "Enable lots of compile warnings" on)

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
	set(BITS32 false)
else()
	set(BITS32 true)
endif()
#option(BITS32 "Build a 32-bit build" off)

option(SHARED_MEM "Enable shared memory" on)

##############################################################################
#                                                                            #
# Compiler specific configuration                                            #
#                                                                            #
##############################################################################

# If this was cmake 2.8.12 or later, then we could use add_compile_options.
# But we can't. The older syntax is really verbose, so use a variable to
# at least make it a bit more compact.
set(CXX_COMPILE_FLAGS "")
set(FORTRAN_COMPILE_FLAGS "")
set(ALL_COMPILE_FLAGS "")
set(ALL_LINKER_FLAGS "")
set(CLUSTER_COMPILE_FLAGS "")
set(CLUSTER_LINK_FLAGS "")

message("Configuring Fortran compilers")
message("-- Build type: ${CMAKE_BUILD_TYPE}")
message("-- Build type: ${CMAKE_BUILD_TYPE}")
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
	message("-- Using ifort")
	# -i_dynamic is incorrectly added in CMake ifort configuration
	set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")
	add_definitions(-D__IFORT)

	set (CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -check bounds")
	set (CLUSTER_COMPILE_FLAGS "-ipo")
	set (CLUSTER_LINK_FLAGS "-ipo")
	
	# Test if we are on Hydra. Because this is an inhomogeneous cluster, it
	# tweaks the available compile options slightly
	site_name(HOSTNAME)
	if (${HOSTNAME} MATCHES "hydra[0-9]*")
		set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -xavx")
	else()
		set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -xHost")
	endif()

elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

	message("-- Using gfortran")
	add_definitions(-D__GFORTRAN__)

	set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -fbounds-check")
	set(FORTRAN_COMPILE_FLAGS "${FORTRAN_COMPILE_FLAGS} -ffree-line-length-none")
	if (NOT BITS32)
		set(FORTRAN_COMPILE_FLAGS "${FORTRAN_COMPILE_FLAGS} -fdefault-real-8 -fdefault-integer-8")
	endif()
	if (WARNINGS)
		set(ALL_COMPILE_FLAGS "${ALL_COMPILE_FLAGS} -Wall -Wextra")
	endif()

else (CMake_Fortran_COMPILER_ID STREQUAL "PGI")

	message("-- Using PGI compilers")
	add_definitions(-D__PGI)

	# Note that pgcpp really is for pre-processing. Not sure why cmake selects
	# this as the c++ compiler - but it doesn't work!
	get_filename_component(CXX_NM ${CMAKE_CXX_COMPILER} NAME)
	if (CXX_NM STREQUAL "pgcpp")
		message("-- Forcing the use of pgc++ not pgcpp")
		set(CMAKE_CXX_COMPILER pgc++)
	endif()

	set(FORTRAN_COMPILE_FLAGS "${FORTRAN_COMPILE_FLAGS} -Msignextend -r8 -i8")
	set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Mbounds")
	set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -tp x64 -fastsse")
	set(ALL_COMPILE_FLAGS "${ALL_COMPILE_FLAGS} -mcmodel=medium")
	set(ALL_LINKER_FLAGS "${ALL_LINKER_FLAGS} -mcmodel=medium")
	if (WARNINGS)
		set(ALL_COMPILE_FLAGS "${ALL_COMPILE_FLAGS} -Minform=warn")
	endif()

	# The linker command -lpthread doesn't work for the PGI linker (due to some
	# internal trickery with the pthread library
	# --> Force the rest of the code to think that libpthreads has already been
	#     found, and then add the magic -pthread command
	# SDS: This no longer seems needed for pgfortran, which does the work automagically
	# set(ALL_LINKER_FLAGS "${ALL_LINKER_FLAGS} -pthread")

else()

	message (FATAL_ERROR "Fortran compiler unknown: ${CMAKE_Fortran_COMPILER_ID}")

endif()

# And options that apply on all compilers
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	add_definitions(-D__DEBUG)
endif()

# Ensure that we signal to the code and the compilers that the correct number
# of bits should be used.
if (BITS32)
	if (CMake_Fortran_COMPILER_ID STREQUAL "PGI")
		set(ALL_COMPILE_FLAGS "${ALL_COMPILE_FLAGS} -pc=32")
		set(ALL_LINKER_FLAGS "${ALL_LINKER_FLAGS} -pc=32")
	else()
		set(ALL_COMPILE_FLAGS "${ALL_COMPILE_FLAGS} -m32")
		set(ALL_LINKER_FLAGS "${ALL_LINKER_FLAGS} -m32")
	endif()
else()
	if (CMake_Fortran_COMPILER_ID STREQUAL "PGI")
		set(ALL_COMPILE_FLAGS "${ALL_COMPILE_FLAGS} -pc=64")
		set(ALL_LINKER_FLAGS "${ALL_LINKER_FLAGS} -pc=64")
	else()
		set(ALL_COMPILE_FLAGS "${ALL_COMPILE_FLAGS} -m64")
		set(ALL_LINKER_FLAGS "${ALL_LINKER_FLAGS} -m64")
	endif()
endif()

# And apply the accumulated flags here
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${ALL_COMPILE_FLAGS} ${FORTRAN_COMPILE_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${ALL_COMPILE_FLAGS} ${CXX_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ALL_COMPILE_FLAGS} ${CXX_COMPILE_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ALL_LINKER_FLAGS}")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${ALL_LINKER_FLAGS}")

# Create the "Cluster" build type. This is based on the Release build type, but
# also includes very slow optimisations (such as -ipo) that one wouldn't want
# to do during a development cycle
set(CMAKE_C_Flags_CLUSTER "${CMAKE_C_FLAGS_RELEASE} ${CLUSTER_COMPILE_FLAGS}")
set(CMAKE_CXX_Flags_CLUSTER "${CMAKE_C_FLAGS_RELEASE} ${CLUSTER_COMPILE_FLAGS}")
set(CMAKE_Fortran_FLAGS_CLUSTER "${CMAKE_Fortran_FLAGS_RELEASE} ${CLUSTER_COMPILE_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS_EXE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} ${CLUSTER_LINK_FLAGS}")

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
target_compile_definitions(neci PUBLIC -DHElement_t=real)
target_compile_definitions(dneci PUBLIC -DHElement_t=real -D__DOUBLERUN)
target_compile_definitions(mneci PUBLIC -DHElement_t=real -D__PROG_NUMRUNS)
target_compile_definitions(kneci PUBLIC -DHElement_t=complex -D__CMPLX)
target_compile_definitions(neci_test PUBLIC -DHElement_t=real)
target_compile_definitions(dneci_test PUBLIC -DHElement_t=real -D__DOUBLERUN)
target_compile_definitions(mneci_test PUBLIC -DHElement_t=real -D__PROG_NUMRUNS)
target_compile_definitions(kneci_test PUBLIC -DHElement_t=complex -D__CMPLX)

# We don't want to compile the test executables unless asked for
set_target_properties(neci_test PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(kneci_test PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(mneci_test PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(dneci_test PROPERTIES EXCLUDE_FROM_ALL TRUE)

#
# If we can get the repository details, then do so.
message("VCS configuration")
find_package(Git)
if (${Git_FOUND})
	message("-- Git has been found")
	Git_WC_INFO(${CMAKE_SOURCE_DIR} "DIR")
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
add_definitions(-DHAVE_SSE2)

add_definitions(-DMAXMEM=99999)
add_definitions(-DDSFMT_MEXP=19937)
add_definitions(-D_CONFIG="cmake-${CMAKE_Fortran_COMPILER_ID}-${CMAKE_BUILD_TYPE}")

# Enable shared memory, unless disabled
if (SHARED_MEM)
	add_definitions(-D__SHARED_MEM)
endif()

if (NOT BITS32)
	add_definitions(-D__INT64)
	add_definitions(-DPOINTER8)
endif()

include_directories(BEFORE SYSTEM ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/unit_tests)

#
# Some form of lapack is always needed (be that acml, mkl, LAPACK, ...)
message("Configuring BLAS/LAPACK")

find_package(MKL)
if (${MKL_FOUND})
	message("-- Found MKL")
	set(NECI_LIBRARIES ${MKL_LIBRARIES})
	include_directories(${MKL_INCLUDE_DIR})
	message("-- Using libraries: ${MKL_LIBRARIES}")
else()
	# Don't use the required flag on find_package here. Do the check manually.
	# If we use the required flag, then the finder does a strict lookup for
	# the pthread library - which will fail for PGI builds due to incompatibilities
	# inside the PGI linker. These are worked around above.
	find_package(LAPACK)
	if (${LAPACK_FOUND})
		set(NECI_LIBRARIES ${LAPACK_LIBRARIES})
		message("-- Found BLAS/Lapack")
		message("-- Using libraries: ${LAPACK_LIBRARIES}")
	else()
		message(FATAL_ERROR "-- BLAS/Lapack not found")
	endif()
endif()

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
if (PARALLEL)

	find_package(MPI REQUIRED)
	if (${MPI_FOUND})
		message("-- Using MPI")
		add_definitions(-DPARALLEL)

		# If we are using ifort, and we have not found the mpiifort wrapper, then
		# we need to force things to use the mpiifort wrapper executable, instead of
		# doing the normal stuff.
		#
		# Otherwise FindMPI can end up finding the GNU stuff, e.g. when using
		# mpi.ibm, which has both sets of wrappers avaialable
		get_filename_component(WRAPPER_NM ${MPI_Fortran_COMPILER} NAME)
		if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" AND NOT WRAPPER_NM STREQUAL "mpiifort")
			message("-- Forcing use of mpiifort wrapper")
			set(MPI_Fortran_COMPILER mpiifort)
			set(MPI_C_COMPILER mpiicc)
			set(MPI_CXX_COMPILER mpiicpc)
			find_package(MPI REQUIRED)
			if (NOT ${MPI_FOUND})
				message(FATAL_ERROR "No suitable MPI wrapper found for ifort")
			endif()
			message("-- Fortran compiler: ${MPI_Fortran_COMPILER}")
			message("-- C compiler: ${MPI_C_COMPILER}")
			message("-- C++ compiler: ${MPI_CXX_COMPILER}")
			set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
			set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
			set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
		else()
			set(NECI_LIBRARIES ${NECI_LIBRARIES} ${MPI_Fortran_LIBRARIES})
			set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fotran_COMPILE_FLAGS}")
			include_directories(${MPI_INCLUDE_PATH})
			include_directories(${MPI_Fortran_INCLUDE_PATH})
			message("-- MPI flags (fortran): ${MPI_Fortran_COMPILE_FLAGS}")
			message("-- MPI include: ${MPI_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH}")
			#message("-- MPI link: ${MPI_Fortran_LINK_FLAG}")
			message("-- MPI libraries: ${MPI_Fortran_LIBRARIES}")
		endif()

		# Prevent poisoning our testing system on later runs
		unset(MPI_Fortran_COMPILER)
		unset(MPI_Fortran_COMPILER CACHE)
		unset(MPI_C_COMPILER)
		unset(MPI_C_COMPILER CACHE)
		unset(MPI_CXX_COMPILER)
		unset(MPI_CXX_COMPILER CACHE)
	endif()

	message("Configuring parallel build")
	if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
	endif()
else()
	message("Parallel build disabled")
endif()

#
# HDF5 as an external dependency
if (HDF)
	message("Configuring HDF5 popsfiles")

	# Note that detection of intel MPI fails, as the names of the wrappers fail
	# --> Need to wrap this for HDF compilation
	if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
		message("Tweaking environment vars")
		set(CONFIGURE_OVERRIDE CC=mpiicc FC=mpiifort F9X=mpiifort CXX=mpiicpc)
	elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
		set(CONFIGURE_OVERRIDE CC=mpicc FC=mpif90 F9X=mpif90 CXX=mpic++ CPP=cpp CXXPP=cpp)
	else()
		set(CONFIGURE_OVERRIDE "")
	endif()

	# This allows us to build hdf5, if we so choose.
	set(HDF_DIR ${CMAKE_CURRENT_BINARY_DIR}/hdf5)
	include(ExternalProject)
	ExternalProject_Add(
		hdf5
		# -- Download step ---
		PREFIX ${HDF_DIR}-prefix
		URL https://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.15-patch1/src/hdf5-1.8.15-patch1.tar.gz
		URL_MD5 4467c25ed9c0b126b194a4d9d66c29ac
		# -- Configure step --
		SOURCE_DIR ${HDF_DIR}
		CONFIGURE_COMMAND env ${CONFIGURE_OVERRIDE} ${HDF_DIR}/configure --enable-parallel --enable-fortran --enable-fortran2003 --prefix=${HDF_DIR}
		# -- Build step ------
		BUILD_COMMAND "" #make && make install
		BUILD_IN_SOURCE 1
		# -- install step ----
		INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/hdf5
		# INSTALL_COMMAND "make install"
	)
	set_target_properties(hdf5 PROPERTIES EXCLUDE_FROM_ALL TRUE)
	set(ENV{HDF5_ROOT} ${HDF_DIR})

	# Because we may be changing the location of the HDF5 libraries after building
	# them, we need to clear existing settings for it
	unset(HDF5_hdf5_LIBRARY_RELEASE)
	unset(HDF5_hdf5_LIBRARY_RELEASE CACHE)
	unset(HDF5_hdf5_fortran_LIBRARY_RELEASE)
	unset(HDF5_hdf5_fortran_LIBRARY_RELEASE CACHE)

	#
	# Find an appropriate HDF5 package
	# n.b. the default HDF5 searcher does not check that the fortran module
	#      was produced with a compatible compiler. As such, we test this
	#      manually, and explicitly, by building a test file contained in tools/
	find_package(HDF5 COMPONENTS Fortran)
	if (${HDF5_FOUND})
		execute_process(
			COMMAND	${CMAKE_Fortran_COMPILER} -I ${HDF5_INCLUDE_DIRS} -c ${CMAKE_SOURCE_DIR}/tools/hdf_module_test.f90
			WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
			RESULT_VARIABLE TEST_RES
		)
		if (NOT ${TEST_RES} STREQUAL "0")
			set (HDF5_FOUND false)
			message("-- HDF5 fortran module not compatible")
		endif()

		if (NOT ${HDF5_IS_PARALLEL})
			message("-- HDF5 not built with MPI support")
			set (HDF5_FOUND false)
		endif()
	endif()

	#
	# If we still think the HDF5 module exits, after testing that it functions
	# correctly, then go ahead and use it!
	if (${HDF5_FOUND})
		message("-- Found HDF5")
		message("-- definitions: ${HDF5_DEFINITIONS}")
		message("-- libs: ${HDF5_LIBRARIES} ${HDF5_Fortran_LIBRARIES}")
		message("-- include: ${HDF5_INCLUDE_DIRS}")
		set(NECI_LIBRARIES ${NECI_LIBRARIES} ${HDF5_Fortran_LIBRARIES} ${HDF5_LIBRARIES})
		add_definitions(${HDF5_DEFINITIONS})
		add_definitions(-D__USE_HDF)
		include_directories(${HDF5_INCLUDE_DIRS})
	else()
		message("-- HDF5 not found")
		message("-- To build using HDF5, first run \"make hdf5\" and then rerun cmake")
	endif()

else()
	message("Disabling HDF5 popsfiles")
endif()


# Ensure that all of our executables link correctly
foreach(target ${TARGETS})
	target_link_libraries(${target} ${NECI_LIBRARIES})
endforeach()
