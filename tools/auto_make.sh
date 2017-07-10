#!/bin/sh

# this script calls the correct make targets on the test machine 
# Usage: 
#	auto_make.sh <config_file_name>

echo "Calling make for $@"

if [ "fkf-ifort hdf5" == "$@" ] || [ "gfortran-hdf5" == "$@" ]; then 
	make VERBOSE=1 hdf5
else
	make VERBOSE=1 -j2 neci mneci kneci kmneci dneci test_neci test_mneci test_kneci test_kmneci test_dneci
fi

