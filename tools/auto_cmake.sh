#!/bin/sh

# This script automatically calls cmake with the correct arguments on the test machine
#
# Usage:
#    auto_module.sh <config_file_name>

echo "Calling cmake for: $@"

if [ "gfortran" == "$@" ]; then
	cmake -DENABLE_HDF5=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx ..
elif [ "gfortran-debug" == "$@" ]; then
	cmake -DENABLE_HDF5=ON -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx ..
elif [ "gfortran-debug-nohdf5" == "$@" ]; then
	cmake -DENABLE_HDF5=OFF -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx ..
elif [ "gfortran-doc" == "$@" ]; then
	cmake -DENABLE_HDF5=ON -DENABLE_DOC=ON -DSHORT_FORD_COMPILATION=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx ..
elif [ "gfortran-debug-integer8" == "$@" ]; then
	cmake -DENABLE_HDF5=ON -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_FLAGS="-fdefault-integer-8" -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx ..
elif [ "ifort" == "$@" ]; then
	cmake -DENABLE_HDF5=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=mpiifort -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc ..
elif [ "ifort-debug" == "$@" ]; then
	cmake -DENABLE_HDF5=ON -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=mpiifort -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc ..
elif [ "ifort18" == "$@" ]; then
	cmake -DENABLE_HDF5=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=mpiifort -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc ..
elif [ "ifort18-debug" == "$@" ]; then
	cmake -DENABLE_HDF5=ON -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=mpiifort -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc ..
elif [ "gfortran-self_build_hdf5" == "$@" ]; then
	cmake -DENABLE_BUILD_HDF5=ON -DENABLE_HDF5=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=mpifort -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx ..
else
	echo "Module not executed"
fi
