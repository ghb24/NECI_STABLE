#!/bin/sh

# This script automatically calls cmake with the correct arguments on the test machine
# 
# Usage:
#    auto_module.sh <config_file_name>

echo "Calling cmake for: $@"

if [ "gfortran-simple" == "$@" ] || [ "fkf-ifort" == "$@" ]; then
	cmake -DCMAKE_BUILD_TYPE=Release ..
elif [ "gfortran-simple -g" == "$@" ] || [ "fkf-ifort -g" == "$@" ]; then
	cmake -DCMAKE_BUILD_TYPE=Debug ..
elif [ "pgi-simple" == "$@" ]; then
	cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_SHARED_MEMORY=off ..
elif [ "pgi-simple -g" == "$@" ]; then
	cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_SHARED_MEMORY=off ..
else
	echo "Module not executed"
fi

