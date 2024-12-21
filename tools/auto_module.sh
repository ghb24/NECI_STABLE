#!/bin/sh

# This script provides the module configurations for the test machine on altest
#
# Usage:
#    auto_module.sh <config_file_name>

source /usr/share/Modules/3.2.10/init/sh
export MODULEPATH="${MODULEPATH}:/home/cicdadmin/modules"
module purge

echo "Loading modules for: $@"

if [ "ifort21-debug" == "$@" ] || [ "ifort21" == "$@" ]; then
    export HDF5_ROOT=/opt/hdf-1.12.2_ifort21_mpi
    module load compiler/2023.0.0 mpi/2021.8.0
elif [ "ifort18-debug" == "$@" ] || [ "ifort18" == "$@" ]; then
    export HDF5_ROOT=/opt/hdf-1.12.2_ifort18_mpi
    source /usr/local/server/IntelStudio_2018/parallel_studio_xe_2018.4.057/psxevars.sh intel64
elif [ "gfortran-debug" == "$@" ] || [ "gfortran-fastdebug" == "$@" ] || [ "gfortran" == "$@" ] || [ "gfortran-short-doc" == "$@" ] || [ "gfortran-doc" == "$@" ] || [ "gfortran-debug-integer8" == "$@" ]; then
    module load hdf5/1.14.5
    module load openmpi/3.1.6
elif [ "gfortran-self_build_hdf5" == "$@" ] || [ "gfortran-debug-nohdf5" == "$@" ]; then
    module load openmpi/3.1.6
else
	echo "Module configuration not set"
fi
