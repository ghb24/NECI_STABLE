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

neci_add_option(
    FEATURE BUILD_HDF5
    DEFAULT off
    DESCRIPTION "Build HDF5 in the source tree, for compatibility with the build configuration" )

if ( HAVE_BUILD_HDF5 )

    # If we are searching for HDF5, we want to create the pseudo-target to build hdf5.

	# Note that detection of intel MPI fails, as the names of the wrappers fail
	# --> Need to wrap this for HDF compilation

    # TODO: Use MPI_Fortran_COMPILER and MPI_C_COMPILER

    set( _c_override "CC=${CMAKE_C_COMPILER}" )
    set( _cxx_override "CXX=${CMAKE_CXX_COMPILER}" )
    set( _fort_override "FC=${CMAKE_Fortran_COMPILER}" "F9X=${CMAKE_Fortran_COMPILER}" )

    #elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
    #	set(_configure_override CC=mpicc FC=mpif90 F9X=mpif90 CXX=mpic++ CPP=cpp CXXPP=cpp)
    #elseif (ARCHER_OVERRIDES)
    #	set(_configure_override CC=cc FC=ftn F9X=ftn CXX=CC CPP=cpp)
    #endif()

    # If we have explicit MPI overrides

    if( MPI_Fortran_COMPILER )
        set( _fort_override "FC=${MPI_Fortran_COMPILER}" "F9X=${MPI_Fortran_COMPILER}" )
    endif()
    if( MPI_C_COMPILER )
        set( _c_override "CC=${MPI_C_COMPILER}" )
    endif()
    if( MPI_CXX_COMPILER )
        set( _cxx_override "CXX=${MPI_CXX_COMPILER}" )
    endif()
    set( _configure_override ${_c_override} ${_cxx_override} ${_fort_override})

    message( STATUS "HDF5 compilation overrides: ${_configure_override}" )

    # Dependencies

    find_package( ZLIB REQUIRED )

    # Create the external build projcet

	set(HDF_DIR "${CMAKE_CURRENT_BINARY_DIR}/hdf5")
    # Some systems use ${HDF_DIR}/lib, others use ${HDF_DIR}/lib64.
    # Make sure to consistently use one of them.
    set(HDF_LIB_DIR "lib")
	include(ExternalProject)
	ExternalProject_Add(
		hdf5
		# -- Download step ---
		PREFIX ${HDF_DIR}-prefix
		#URL https://support.hdfgroup.org/ftp/HDF5/current18/src/hdf5-1.8.20.tar.gz
		URL https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.20/src/hdf5-1.8.20.tar.gz
		URL_MD5 7f2d3fd67106968eb45d133f5a22150f

	# previous 1.8.19 hash:
		# URL https://support.hdfgroup.org/ftp/HDF5/current18/src/hdf5-1.8.19.tar.gz
		# URL_MD5 7f568e2464d4ab0a74d16b23956d900b
	# previous 1.8.18 hash:
	# URL_MD5 dd2148b740713ca0295442ec683d7b1c
        # URL http://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.15-patch1/src/hdf5-1.8.15-patch1.tar.gz
		# URL_MD5 4467c25ed9c0b126b194a4d9d66c29ac
		# -- Configure step --
		SOURCE_DIR ${HDF_DIR}
		CONFIGURE_COMMAND env ${_configure_override} ${HDF_DIR}/configure --enable-parallel --enable-fortran --enable-fortran2003 --prefix=${HDF_DIR} --libdir=${HDF_DIR}/${HDF_LIB_DIR}
		# -- Build step ------
		BUILD_COMMAND "" #make && make install
		BUILD_IN_SOURCE 1
		# -- install step ----
		INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/hdf5
		# INSTALL_COMMAND "make install"
	)

    # Add this target to the dependencies of evrything else, so it is automagically built first

    message( STATUS "Adding global dependency to hdf5 target" )
    set( ${PROJECT_NAME}_GLOBAL_DEPENDENCIES ${${PROJECT_NAME}_GLOBAL_DEPENDENCIES} hdf5 )

    # Add the appropriate variable components

    set( HDF5_FOUND on )
    set( HDF5_LIBRARIES ${HDF_DIR}/${HDF_LIB_DIR}/libhdf5.a dl ${ZLIB_LIBRARIES} )
    set( HDF5_Fortran_LIBRARIES ${HDF_DIR}/${HDF_LIB_DIR}/libhdf5_fortran.a )
    set( HDF5_INCLUDE_DIRS ${HDF_DIR}/include )
    set( HDF5_Fortran_INCLUDE_DIRS )

else() # Not building hdf5 ...

	#
	# Find an appropriate HDF5 package
	# n.b. the default HDF5 searcher does not check that the fortran module
	#      was produced with a compatible compiler. As such, we test this
	#      manually, and explicitly, by building a test file contained in tools/

    # If NECI is a git submodule of other programs e.g. MOLCAS,
    # that don't use the Fortran compoments, we want to make sure,
    # that find_package does not use the cached values, but really makes
    # a new search.
    set (HDF5_FOUND OFF)
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
endif()


if ( NOT HDF5_NECI_FOUND )

    set( HDF5_NECI_FOUND ${HDF5_FOUND} )
    if ( HDF5_FOUND )
        set( HDF5_NECI_INCLUDE_PATH ${HDF5_Fortran_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
        set( HDF5_NECI_DEFINITIONS ${HDF5_DEFINITIONS})
        if (NOT HDF5_NECI_FIND_QUIETLY)
            message(STATUS "HDF5 found")
        endif()

		# Filter for debug/optimized bits in libraries
		# If this filtering becomes useful elsewhere, it should be extracted into a macro.
		set( HDF5_NECI_LIBRARIES "" )
		set( _skipnext OFF )
		foreach( _lib ${HDF5_Fortran_LIBRARIES} ${HDF5_LIBRARIES} )
			if( _lib STREQUAL "debug" )
				if( NOT CMAKE_BUILD_TYPE STREQUAL "DEBUG" )
					set( _skipnext ON )
				endif()
			elseif( _lib STREQUAL "optimized" )
				if( CMAKE_BUILD_TYPE STREQUAL "DEBUG" )
					set( _skipnext ON )
				endif()
			else()
				if( _skipnext )
					set( _skipnext OFF )
				else()
					list( APPEND HDF5_NECI_LIBRARIES ${_lib} )
				endif()
			endif()
		endforeach()

    else()
        if ( HDF5_NECI_FIND_REQUIRED )
            message(FATAL_ERROR "Package HDF5 required but not found")
        else()
            if (NOT HDF5_NECI_FIND_QUIETLY)
                message(STATUS "Package HDF5 not found")
            endif()
        endif()
    endif()

endif()

