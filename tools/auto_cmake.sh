#!/bin/sh

# This script automatically calls cmake with the correct arguments on the test machine
#
# Usage:
#    auto_module.sh <config_file_name>

BUILD_TYPE=$1
SOURCE_DIR=${2:-..}

echo "Calling cmake for: $BUILD_TYPE"

if [ "gfort3" == "$BUILD_TYPE" ] || [ "gfort5" == "$BUILD_TYPE" ]; then
	cmake -DENABLE_HDF5=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx  $SOURCE_DIR
elif [ "gfort3-debug" == "$BUILD_TYPE" ] || [ "gfort5-debug" == "$BUILD_TYPE" ]; then
	cmake -DENABLE_HDF5=ON -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx  $SOURCE_DIR
elif [ "gfort3-fastdebug" == "$BUILD_TYPE" ] || [ "gfort5-fastdebug" == "$BUILD_TYPE" ]; then
	cmake -DENABLE_HDF5=ON -DCMAKE_BUILD_TYPE=FastDebug -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx  $SOURCE_DIR
elif [ "gfort3-debug-nohdf5" == "$BUILD_TYPE" ] || [ "gfort5-debug-nohdf5" == "$BUILD_TYPE" ]; then
	cmake -DENABLE_HDF5=OFF -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx  $SOURCE_DIR
elif [ "gfort3-doc" == "$BUILD_TYPE" ] || [ "gfort5-doc" == "$BUILD_TYPE" ]; then
	cmake -DENABLE_HDF5=ON -DENABLE_DOC=ON -DSHORT_FORD_COMPILATION=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx  $SOURCE_DIR
elif [ "gfort3-short-doc" == "$BUILD_TYPE" ] || [ "gfort5-short-doc" == "$BUILD_TYPE" ]; then
	cmake -DENABLE_HDF5=ON -DENABLE_DOC=ON -DSHORT_FORD_COMPILATION=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx  $SOURCE_DIR
elif [ "gfort3-debug-integer8" == "$BUILD_TYPE" ] || [ "gfort5-debug-integer8" == "$BUILD_TYPE" ]; then
	cmake -DENABLE_HDF5=ON -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_FLAGS="-fdefault-integer-8" -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx  $SOURCE_DIR
elif [ "ifort" == "$BUILD_TYPE" ]; then
	cmake -DENABLE_HDF5=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=mpiifort -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc  $SOURCE_DIR
elif [ "ifort-debug" == "$BUILD_TYPE" ]; then
	cmake -DENABLE_HDF5=ON -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=mpiifort -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc  $SOURCE_DIR
elif [ "ifort18" == "$BUILD_TYPE" ] || [ "ifort21" == "$BUILD_TYPE" ]; then
	cmake -DENABLE_HDF5=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=mpiifort -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc  $SOURCE_DIR
elif [ "ifort18-debug" == "$BUILD_TYPE" ] || [ "ifort21-debug" == "$BUILD_TYPE" ]; then
	cmake -DENABLE_HDF5=ON -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=mpiifort -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc  $SOURCE_DIR
elif [ "gfort3-self_build_hdf5" == "$BUILD_TYPE" ] || [ "gfort5-self_build_hdf5" == "$BUILD_TYPE" ]; then
	cmake -DENABLE_BUILD_HDF5=ON -DENABLE_HDF5=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx  $SOURCE_DIR
else
	echo "Module not executed"
fi
