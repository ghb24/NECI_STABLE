#!/bin/sh

# This script provides the module configurations for the test machine on altest
# 
# Usage:
#    auto_module.sh <config_file_name>

echo "Loading modules for: $@"

if [ "fkf-ifort -g" == "$@" ] || [ "fkf-ifort" == "$@" ] || [ "fkf-ifort-new" == "$@" ]; then
        export HDF5_ROOT=/usr/lib/custom_hdf5_ifort
        module load ifort/18.0.1 mpi.intel/5.0.3 #hdf5-intel
elif [ "gfortran-simple -g" == "$@" ] || [ "gfortran-simple" == "$@" ]; then
        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64/mpi/gcc/openmpi3/lib64
	export HDF5_ROOT=/usr/lib/custom_hdf5_gfortran
#	module load gnu-openmpi/3.0.0 hdf5-gfortran/1.8.20
elif [ "pgi-simple -g" == "$@" ] || [ "pgi-simple" == "$@" ]; then
	module load PrgEnv-pgi/15.4 hdf5-pgi
elif [ "fkf-ifort-hdf5" == "$@" ]; then 
	module load ifort/18.0.1 mpi.intel/5.0.3
elif [ "gfortran-hdf5" == "$@" ]; then
        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64/mpi/gcc/openmpi3/lib64
#	module load gnu-openmpi/3.0.0
elif [ "fkf-ifort-latest" == "$@" ]; then 
	module load fkf-ifort mpi.intel
else
	echo "Module configuration not set"
fi

