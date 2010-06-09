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

int main ()
{
	// Get the current working directory and substitute as appropriate
	char temp[MAXPATHLEN];
	string cwd;
	if (getcwd(temp, MAXPATHLEN)) {
		cwd = string(temp);
		std::replace (cwd.begin(), cwd.end(), '/', '_');
	} else {
		cerr << "Error getting current path\n";
		return 1;
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
	cout << "Searching for residual shared memory devices... " << endl;
	while ((dirp = readdir(dp)) != NULL) {
		string nm (dirp->d_name);
		// Test if this is a shared memory device originating from an
		// instantiation of NECI in the current directory, and if so that
		// it is not one from a subdirectory (which may still be legitimately
		// running)
		if (nm.length() > cwd.length() &&
			nm.substr(0, cwd.length()).compare(cwd) == 0 &&
			nm[cwd.length()] != '_') {

			cout << "Removing: " << nm << endl;
			shm_unlink (nm.c_str());
		}
	}

	closedir (dp);
	return 0;
}
