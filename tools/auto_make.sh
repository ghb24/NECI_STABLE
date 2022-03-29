#!/bin/sh

# this script calls the correct make targets on the test machine
# Usage:
#	auto_make.sh <config_file_name>

echo "Calling make for $@"

if [ "gfortran-self_build_hdf5" == "$@" ]; then
	make VERBOSE=1 hdf5

	# very ugle workaround to get this compilation working on altest..
	# for some reason on altest the hdf5 libs get installed in hdf5/lib64 instead of hdf5/lib
	# have to solve this more elegantly
	if [ -d "hdf5/lib" ]; then
		make -j2 neci
	elif [ -d "hdf5/lib64" ]; then
		mkdir hdf5/lib
		cp hdf5/lib64/* hdf5/lib
		make -j2 neci mneci kneci kmneci dneci test_neci test_mneci test_kneci test_kmneci test_dneci
	fi

else
	make VERBOSE=1 -j2 neci mneci kneci kmneci dneci test_neci test_mneci test_kneci test_kmneci test_dneci
fi

