#!/bin/sh

# This script provides the module configurations for the test machine on altest
#
# Usage:
#    auto_module.sh <config_file_name>

echo "Loading modules for: $@"

if [ "ifort-debug" == "$@" ] || [ "ifort" == "$@" ]; then
    export HDF5_ROOT=/opt/hdf-1.8.20_ifort_19
    module load ifort/19.1.1 mpi.intel/2019.7
elif [ "gfortran-debug" == "$@" ] || [ "gfortran" == "$@" ] || [ "gfortran-integer8" == "$@" ]; then
    export HDF5_ROOT=/opt/hdf-1.8.20_gfort_19
    module load gnu-openmpi/3.1.6
elif [ "gfortran-self_build_hdf5" == "$@" ]; then
    module load gnu-openmpi/3.1.6
else
	echo "Module configuration not set"
fi

