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
elif [ "gfort7-debug" == "$@" ] || [ "gfort7-fastdebug" == "$@" ] || [ "gfort7" == "$@" ] || [ "gfort7-short-doc" == "$@" ] || [ "gfort7-doc" == "$@" ] || [ "gfort7-debug-integer8" == "$@" ]; then
    module load hdf5/1.14.5
    module load openmpi/3.1.6
elif [ "gfort14-debug" == "$@" ] || [ "gfort14-fastdebug" == "$@" ] || [ "gfort14" == "$@" ] || [ "gfort14-short-doc" == "$@" ] || [ "gfort14-doc" == "$@" ] || [ "gfort14-debug-integer8" == "$@" ]; then
    module load hdf5/1.14.5
    module load openmpi/5.0.6
elif [ "gfort7-self_build_hdf5" == "$@" ] || [ "gfort7-debug-nohdf5" == "$@" ]; then
    module load openmpi/3.1.6
elif [ "gfort14-self_build_hdf5" == "$@" ] || [ "gfort14-debug-nohdf5" == "$@" ]; then
    module load openmpi/5.0.6
else
	echo "Module configuration not set"
fi
