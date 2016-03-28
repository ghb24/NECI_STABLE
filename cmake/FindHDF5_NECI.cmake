# Instructs CMake to find an appropriate HDF5 package (or build it if required)
#
# i) CMake will examine the path set in the environment variable HDF5_ROOT (which is set by the module
#    system in important locations).
#
# ii) This creates the target "hdf5" which will download and build a version of hdf5 in an hdf5/ subdir
#     of the build directory. This hdf5 module can then be found by running cmake again.
#
# Provides:
# HDF5_NECI_FOUND
# HDF5_NECI_INCLUDE_PATHS
# HDF5_NECI_LIBRARIES
# HDF5_NECI_DEFINITIONS

if ( NOT HDF5_NECI_FOUND )

    # -------------------------------------------------------------------------------------------------------

    # If we are searching for HDF5, we want to create the pseudo-target to build hdf5.

	# Note that detection of intel MPI fails, as the names of the wrappers fail
	# --> Need to wrap this for HDF compilation
	if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
		set(_configure_override CC=mpiicc FC=mpiifort F9X=mpiifort CXX=mpiicpc)
	elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
		set(_configure_override CC=mpicc FC=mpif90 F9X=mpif90 CXX=mpic++ CPP=cpp CXXPP=cpp)
	elseif (ARCHER_OVERRIDES)
		set(_configure_override CC=cc FC=ftn F9X=ftn CXX=CC CPP=cpp)
	else()
		set(_configure_override "")
	endif()

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
	set(SYS_HDF5_ROOT $ENV{HDF5_ROOT})
	set(ENV{HDF5_ROOT} ${HDF_DIR}:$ENV{HDF5_ROOT})
	find_package(HDF5 COMPONENTS Fortran)
	if (${HDF5_FOUND})
		execute_process(
			COMMAND	${CMAKE_Fortran_COMPILER} -I ${HDF5_INCLUDE_DIRS} -c ${PROJECT_SOURCE_DIR}/tools/hdf_module_test.f90
			WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
			RESULT_VARIABLE TEST_RES
		)
		if (NOT ${TEST_RES} STREQUAL "0")
			set (HDF5_FOUND false)
			message(WARNING "HDF5 fortran module not compatible")
		endif()

		if (NOT ${HDF5_IS_PARALLEL})
			message(WARNING "HDF5 not built with MPI support")
			set (HDF5_FOUND false)
		endif()
	endif()
	set(ENV{HDF5_ROOT} ${SYS_HDF5_ROOT})


# TODO: Add special message if HDF5 has not been found
# list(APPEND ${PROJECT_NAME}_SPECIAL_MESSAGES 

    set( HDF5_NECI_FOUND ${HDF5_FOUND} )
    if ( HDF5_FOUND )
        set( HDF5_NECI_LIBRARIES ${HDF5_Fortran_LIBRARIES} ${HDF5_LIBRARIES} )
        set( HDF5_NECI_INCLUDE_PATH ${HDF5_Fortran_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
        set( HDF5_NECI_DEFINITIONS ${HDF5_DEFINITIONS})
        if (NOT HDF5_NECI_FIND_QUIETLY)
            message(STATUS "HDF5 found")
        endif()
    else()
        if ( HDF5_NECI_FIND_REQUIRED )
            message(FATAL_ERROR "Package HDF5 required but not found")
        else()
            if (NOT HDF5_NECI_FIND_QUIETLY)
                message(STATUS "Package HDF5 not found")
            endif()
        endif()
    endif()

##    find_package( MPI )
##
##    # Map the output of the find to the MPI_NECI finder
##    set( MPI_NECI_FOUND ${MPI_FOUND})
##    if ( MPI_FOUND )
##        set(MPI_NECI_LIBRARIES ${MPI_LIBRARIES})
##        set(MPI_NECI_Fortran_LIBRARIES ${MPI_Fortran_LIBRARIES})
##        set(MPI_NECI_INCLUDE_PATH ${MPI_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH})
##        if (NOT MPI_NECI_FIND_QUIETLY )
##            message(STATUS "MPI found")
##        endif()
##    else()
##        if ( MPI_NECI_FIND_REQUIRED )
##            message(FATAL_ERROR "Package MPI required, but not found")
##        else()
##            if (NOT MPI_NECI_FIND_QUIETLY )
##                message(STATUS "Package MPI not found")
##            endif()
##        endif()
##    endif()

endif()

