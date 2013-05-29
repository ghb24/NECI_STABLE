#include <iostream>
#include <string>
#include <algorithm>
#include <sys/types.h>
#include <sys/param.h>
#include <sys/mman.h>
#include <dirent.h>
#include <errno.h>

using std::string;
using std::cerr;
using std::cout;
using std::endl;

// Remove any lingering shared memory devices (e.g. a careless programmer has
// not deallocated memory correctly or caused NECI to crash).
// Note that this *must* be run from the same working directory as the
// instantiation of NECI which created the shared memory devices.
int main (int argc, char* argv[])
{
	// If we specify all, then we don't need the directory length etc.
	bool del_all = false;
	if (argc > 1 && string(argv[1]) == "all") {
		cout << "Removing all NECI-like shared memory devices\n";
		del_all = true;
	}

	// Get the current working directory and substitute as appropriate
	char temp[MAXPATHLEN];
	string cwd;
	if (!del_all) {
		if (getcwd(temp, MAXPATHLEN)) {
			cwd = string(temp);
			std::replace (cwd.begin(), cwd.end(), '/', '_');
		} else {
			cerr << "Error getting current path\n";
			return 1;
		}
	}

	// Open the shared memory directory
	DIR * dp;
	if ((dp = opendir ("/dev/shm")) == NULL) {
		cerr << "Error opening shared memory device directory /dev/shm\n";
		return errno;
	}

	// Iterate through the devices, and see if there are any we are
	// interested in
	struct dirent * dirp;
	bool found_shm = false;
	cout << "Searching for residual shared memory devices... " << endl;
	while ((dirp = readdir(dp)) != NULL) {
		string nm (dirp->d_name);
		// Test if this is a shared memory device originating from an
		// instantiation of NECI in the current directory, and if so that
		// it is not one from a subdirectory (which may still be legitimately
		// running)
		// NECI creates shared memory devices with the desired name (e.g.
		// "umat") prepended by the working directory with / replaced by _.
		if (nm[0] == '_') {
			if (del_all ||
				(nm.length() > cwd.length() &&
				 nm.substr(0, cwd.length()).compare(cwd) == 0 &&
				 nm[cwd.length()] != '_')) {

				cout << "Removing: " << nm << endl;
				shm_unlink (nm.c_str());
				found_shm = true;
			}
		}
	}

	if (! found_shm) {
		cout << "No residual shared memory devices found." << endl;
	}

	closedir (dp);
	return 0;
}
