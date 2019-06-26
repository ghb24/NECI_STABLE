#!/bin/sh

# This script automatically calls cmake with the correct arguments on the test machine
# 
# Usage:
#    auto_module.sh <config_file_name>

echo "Calling cmake for: $@"

if [ "gfortran-simple" == "$@" ]; then
	cmake -DENABLE_HDF5=ON -DCMAKE_BUILD_TYPE=Release ..
elif [ "fkf-ifort" == "$@" ] || [ "fkf-ifort-new" == "$@" ]; then
	cmake -DENABLE_HDF5=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=mpiifort -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc ..
elif [ "gfortran-simple -g" == "$@" ]; then
	cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_HDF5=ON ..
elif [ "fkf-ifort -g" == "$@" ]; then
	cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_HDF5=ON -DCMAKE_Fortran_COMPILER=mpiifort -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc ..
elif [ "pgi-simple" == "$@" ]; then
	cmake -DENABLE_HDF5=ON -DCMAKE_BUILD_TYPE=Release -DENABLE_SHARED_MEMORY=off ..
elif [ "pgi-simple -g" == "$@" ]; then
	cmake -DENABLE_HDF5=ON -DCMAKE_BUILD_TYPE=Debug -DENABLE_SHARED_MEMORY=off ..
elif [ "fkf-ifort-hdf5" == "$@" ]; then
	cmake -DENABLE_BUILD_HDF5=ON -DENABLE_HDF5=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=mpiifort -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc ..
elif [ "gfortran-hdf5" == "$@" ]; then 
	cmake -DENABLE_BUILD_HDF5=ON -DENABLE_HDF5=ON -DCMAKE_BUILD_TYPE=Release ..
else
	echo "Module not executed"
fi

