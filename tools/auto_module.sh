#!/bin/sh

# This script provides the module configurations for the test machine on altest
# 
# Usage:
#    auto_module.sh <config_file_name>

echo "Loading modules for: $@"

if [ "$@" == "fkf-ifort" ] || [ "$@" == "fkf-ifort -g" ]; then
	module load ifort/15.0.3 mpi.intel/5.0.3
else
	echo "Module configuration not set"
fi
	
