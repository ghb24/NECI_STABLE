#!/bin/sh

# this script calls the correct make targets on the test machine
# Usage:
#	auto_make.sh <config_file_name>
# cmake --build . --target neci test_neci -j 6 -v --
echo "Calling make for $@"

if [ "gfortran-self_build_hdf5" == "$@" ]; then

	make VERBOSE=1 hdf5

	# very ugly workaround to get this compilation working on altest..
	# for some reason on altest the hdf5 libs get installed in hdf5/lib64 instead of hdf5/lib
	# have to solve this more elegantly
	if [ -d "hdf5/lib" ]; then
        # make -j 2
		cmake --build . -j 2 --
	elif [ -d "hdf5/lib64" ]; then
		mkdir hdf5/lib
		cp hdf5/lib64/* hdf5/lib
		# make -j 2
		cmake --build . -j 2 --
	fi

elif [ "gfortran-doc" == "$@" ] || [ "gfortran-short-doc" == "$@" ]; then
    # make VERBOSE=1 -j 2 doc
    cmake --build . -j 2 --target doc -v --
else
	# make VERBOSE=1 -j 2
    cmake --build . -j 2 -v --
fi
