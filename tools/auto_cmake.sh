#!/bin/sh

# This script automatically calls cmake with the correct arguments on the test machine
# 
# Usage:
#    auto_module.sh <config_file_name>

echo "Calling cmake for: $@"

if [ "pgi-simple" == "$@" ] || [ "gfortran-simple" == "$@" ] || [ "fkf-ifort" == "$@" ]; then
	cmake -DCMAKE_BUILD_TYPE=Release ..
elif [ "pgi-simple -g" == "$@" ] || [ "gfortran-simple -g" == "$@" ] || [ "fkf-ifort -g" == "$@" ]; then
	cmake -DCMAKE_BUILD_TYPE=Debug ..
else
	echo "Module not executed"
fi

