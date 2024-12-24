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
elif [ "gfort3-debug" == "$@" ] || [ "gfort3-fastdebug" == "$@" ] || [ "gfort3" == "$@" ] || [ "gfort3-short-doc" == "$@" ] || [ "gfort3-doc" == "$@" ] || [ "gfort3-debug-integer8" == "$@" ]; then
    module load hdf5/1.14.5
    module load openmpi/3.1.6
elif [ "gfort5-debug" == "$@" ] || [ "gfort5-fastdebug" == "$@" ] || [ "gfort5" == "$@" ] || [ "gfort5-short-doc" == "$@" ] || [ "gfort5-doc" == "$@" ] || [ "gfort5-debug-integer8" == "$@" ]; then
    module load hdf5/1.14.5
    module load openmpi/5.0.6
elif [ "gfort3-self_build_hdf5" == "$@" ] || [ "gfort3-debug-nohdf5" == "$@" ]; then
    module load openmpi/3.1.6
elif [ "gfort5-self_build_hdf5" == "$@" ] || [ "gfort5-debug-nohdf5" == "$@" ]; then
    module load openmpi/5.0.6
else
	echo "Module configuration not set"
fi
