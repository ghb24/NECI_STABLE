#include <stdio.h>
#include <stdint.h>
#include <stdarg.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <fcntl.h>
#include <map>
#include <list>
#include <string>
#include <fstream>
#include <iostream>
#include <stdarg.h>

#ifdef _WIN32
#include <windows.h>
#include <tchar.h>
#else
#include <sys/mman.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/param.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#endif

using std::map;
using std::string;
using std::list;
using std::ofstream;

extern "C" void stop_all_c (const char* a, const char* b);
extern "C" void mpibarrier_c (int * error);
extern "C" void print_cstr (const char* s);

void fort_printf (const char * fmt, ...)
{
	// Create a buffer to process the printf arguments.
	const int buflen = 1024;
	char buf[buflen];

	// Process the arguments
	va_list args;
	va_start(args, fmt);
	vsnprintf (buf, buflen, fmt, args);
	va_end(args);

	// Ensure that the string is terminated, even if it has been truncated
	buf[buflen-1] = '\0';

	// And output via fortran
	print_cstr (buf);
}

// This shared memory mapping is the only bit of this file which is c++
// rather than C, but it is definitely a better way of doing it than 
// writing a map in C!
//
// Store a map with the required details to unlink the shared memory. This
// means that the fortran program only has to keep track of the pointer (which
// it is directly using).
class map_det_t {
public:
	map_det_t () : size(0) {}
#ifdef _WIN32
	map_det_t (string name_, HANDLE hMap_) : name(name_), hMap(hMap_) {}
	HANDLE hMap;
#else
	map_det_t (string name_, int size_) : name(name_), size(size_) {}
#endif
	string name;
	size_t size;
};
map<void*,map_det_t> g_shared_mem_map;
list<void*> g_shm_list;

std::string cwd_name;
std::wstring wcwd_name;

#ifndef _WIN32
//
//
// Test if shared memory is going to be OK. If not, print a warning...
// Only print the warnings on the first run through.
extern "C" bool test_shared_permissions ()
{
	const char* shm_dir = "/dev/shm";
	static bool avail = true;
	static bool tested = false;

	if (!tested) {
		int stat = access (shm_dir, F_OK);
		if (stat) {
			fort_printf ("SHM directory not found: %s\n", strerror(errno));
			avail = false;
		}

		stat = access (shm_dir, R_OK);
		if (stat) {
			fort_printf ("Read access to SHM directory unavailable: %s\n",
					strerror(errno));
			avail = false;
		}

		stat = access (shm_dir, W_OK);
		if (stat) {
			fort_printf ("Write access to SHM directory unavailable: %s\n",
					strerror(errno));
			avail = false;
		}

		// Indicate that the above are warnings, not errors.
		if (!avail)
			fort_printf ("Falling back on System V Shared memory.\n");

		fflush (stdout);
		tested = true;
	}

	return avail;
}
#endif

#ifdef _WIN32

//
// Allocate shared memory in windows
void allocate_shared_windows (const char * name, void ** ptr,
                              const size_t size, int proc)
{
	// Get the current working directory name if needed.
	// --> Prepend to names asked for, to give unique name to avoid collisions.
	// This really ought be done using TCHARs on windows...
	if (wcwd_name.empty()) {
		WCHAR szDir[MAX_PATH] = L"";
		if (!::GetCurrentDirectoryW(MAX_PATH, szDir))
			stop_all_c (__FUNCTION__,
			            "Obtaining current working directory failed");

		wcwd_name = szDir;
		std::replace (wcwd_name.begin(), wcwd_name.end(), L'/', L'_');
		std::replace (wcwd_name.begin(), wcwd_name.end(), L'\\', L'_');
	}

	// Name to use for shared memory.
	WCHAR szNm[MAX_PATH] = L"";
	MultiByteToWideChar(CP_ACP, 0, name, -1, szNm, MAX_PATH);
	std::wstring shared_name = wcwd_name + szNm;

	// Make sure we use the unicode functions. This allows us to deal with
	// more interesting directory names (although not file names directly,
	// as 'name' passed from FORTRAN is somewhat restricted).
	HANDLE hMapFile = NULL;
	if (0 == proc) {
		// One process has to create the mapping
		hMapFile = CreateFileMappingW (INVALID_HANDLE_VALUE, NULL,
		                               PAGE_READWRITE, 0, size,
		                               shared_name.c_str());
		if (NULL == hMapFile)
			stop_all_c (__FUNCTION__,
			            "Could not create file mapping object");
	}

	// All processes need to wait until the object has been created.
	int ierr;
	mpibarrier_c (&ierr);

	if (0 != proc) {
		// The other processes then have to find it.
		hMapFile = OpenFileMappingW (FILE_MAP_ALL_ACCESS, FALSE,
		                             shared_name.c_str());
		if (NULL == hMapFile)
			stop_all_c (__FUNCTION__,
			            "Could not open file mapping object");
	}

	// And again we wait.
	mpibarrier_c (&ierr);

	// Map the object into memory.
	*ptr = MapViewOfFile (hMapFile, FILE_MAP_ALL_ACCESS, 0, 0, size);
	if (NULL == *ptr) {
		CloseHandle(hMapFile);
		stop_all_c (__FUNCTION__,
					"Memory mapping file object failed.");
	}
	
	// Store the handle for later cleanup.
	g_shared_mem_map[*ptr] = map_det_t(name, hMapFile);
	mpibarrier_c (&ierr);
}

#else // _WIN32

//
// Allocate shared memory using POSIX shared memory semantics.
void allocate_shared_posix (const char * name, void ** ptr, const size_t size)
{
	// Get the current working directory name if needed
	// --> Prepend to names asked for, to give unique name to avoid collisions
	if (cwd_name.empty()) {
		char temp[MAXPATHLEN];
		cwd_name = getcwd(temp, MAXPATHLEN) ? string(temp)
			                                : string("NECI_LONG_PATH");
		std::replace (cwd_name.begin(), cwd_name.end(), '/', '_');
	}

	// Acquire a named shared memory file descriptor. nb. Enforce 10 character
	// name limit.
	string shared_name = cwd_name + name;
	int fd = shm_open (shared_name.c_str(), O_CREAT | O_RDWR, (mode_t)0600);
	if (fd == -1)
		stop_all_c (__FUNCTION__,
		          (string("creating shared memory object failed: ") 
		           + strerror(errno)).c_str());

	// Set the length of the named region to the required length
	if (ftruncate(fd, size) == -1)
		stop_all_c (__FUNCTION__,
		          (string("Setting size of shared memory region failed: ") 
		           + strerror(errno)).c_str());

	// Map the region into the current address space.
	*ptr = mmap (NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, 
			fd, 0);
	if (*ptr == MAP_FAILED)
		stop_all_c (__FUNCTION__,
		          (string("Mapping shared memory failed: ")
		           + strerror(errno)).c_str());

	// Once we have mapped the region, it will remain available even when the
	// file descriptor has closed.
	close(fd);

	g_shared_mem_map[*ptr] = map_det_t(shared_name, size);

	// Wait until all threads have reached this point, and then unlink the
	// shared memory object --> Operating system will clean up correctly if
	// the processes bail.
	int ierr;
	mpibarrier_c (&ierr);
	shm_unlink (shared_name.c_str());
}

//
// Allocate shared memory using unix systemV semantics.
void allocate_shared_systemV (const char * name, void ** ptr,
                              const size_t size)
{
	// What path should we use for the shared memory object file?
	string path = string("._NECI_shm$") + string(name);

	// Ensure that the above file exists
	ofstream ftmp (path.c_str());
	ftmp.close();

	// Get the shm key associated with the above file, and the character '~'
	key_t shm_key = ftok (path.c_str(), '~');
	if (shm_key == -1)
		stop_all_c (__FUNCTION__, (string("Error getting System V IPC key: ") 
		                        + strerror(errno)).c_str());

	// Get the shm object reference.
	int shm_id = shmget (shm_key, size, 0644 | IPC_CREAT);
	if (shm_id == -1)
		stop_all_c (__FUNCTION__,
		          (string("Error obtaining shared memory segment: ")
		           + strerror(errno)).c_str());

	// Map shm object into memory
	*ptr = shmat (shm_id, NULL, 0);
	if (*ptr == (void*)-1)
		stop_all_c (__FUNCTION__, (string("Error mapping shared memory: ")
		                         + strerror(errno)).c_str());

	// Wait until all threads have reached this point, and then remove the
	// shared memory control object --> Operating system will clean up
	// correctly if the processes bail.
	int ierr;
	mpibarrier_c (&ierr);
	shmctl (shm_id, IPC_RMID, NULL);
	g_shm_list.push_back(*ptr);
}

#endif // _WIN32

//
//
// This function acquires a named region of shared memory, of the specified
// size, allocating or resizing it if necessary. It then returns this as a
// standard C pointer --> requires fortran 2003 features to use.
extern "C" void alloc_shared_worker (const char * name, void ** ptr,
		                             const size_t size, int proc)
{
#ifdef _WIN32
	allocate_shared_windows (name, ptr, size, proc);
#else
	if (test_shared_permissions())
		allocate_shared_posix (name, ptr, size);
	else
		allocate_shared_systemV (name, ptr, size);
#endif
}

//
// Given a C pointer to a region of shared memory, retrieve all of the
// relevant details from the map, and then deallocate the memory and remove
// the mapping.
extern "C" void dealloc_shared_worker (void * ptr)
{

#ifdef _WIN32
	// Find the shared memroy in the global list (windows)
	map<void*,map_det_t>::iterator det;
	det = g_shared_mem_map.find(ptr);

	if (det == g_shared_mem_map.end())
		stop_all_c (__FUNCTION__,
		            "The specified shared memory was not found.");

	// Unmap and close the shared segment.
	UnmapViewOfFile (ptr);
	CloseHandle(det->second.hMap);

	g_shared_mem_map.erase(det);
#else
	if (test_shared_permissions()) {
		//
		// Deallocate POSIX shared memory
		
		// Find the shared memory in the global list.
		map<void*,map_det_t>::iterator det;
		det = g_shared_mem_map.find(ptr);

		if (det == g_shared_mem_map.end())
			stop_all_c (__FUNCTION__,
			          "The specified shared memory was not found.");

		munmap (ptr, det->second.size);

		g_shared_mem_map.erase(det);
	} else {
		//
		// Deallocate systemV shared memory
		list<void*>::iterator item;
		item = find (g_shm_list.begin(), g_shm_list.end(), ptr);
		if (item == g_shm_list.end())
			stop_all_c (__FUNCTION__,
			          "The specified shared memory was not found.");

		shmdt (ptr);
		g_shm_list.erase(item);
	}
#endif
}


//
// Clean up any shared allocations which have not been properly deallocated.
extern "C" void cleanup_shared_alloc ()
{
#ifdef _WIN32
	// Iterate through the list of shared allocations and clear up (windows)
	map<void*,map_det_t>::iterator iter;
	for (iter = g_shared_mem_map.begin(); iter!= g_shared_mem_map.end();
	     ++iter) {
		size_t size = iter->second.size;
		void * ptr = iter->first;
		string name = iter->second.name;
		fort_printf ("Non-deallocated shared memory found: %s, %d bytes\n", 
		        name.c_str(), int(size));
		UnmapViewOfFile (ptr);
		CloseHandle (iter->second.hMap);
	}
	g_shared_mem_map.clear();

#else
	if (test_shared_permissions()) {
		//
		// Clean up POSIX shared memory

		// Iterate through the list of shared allocations
		map<void*,map_det_t>::iterator iter;
		for (iter = g_shared_mem_map.begin(); iter != g_shared_mem_map.end();
			 ++iter) {
			size_t size = iter->second.size;
			void * ptr = iter->first;
			string name = iter->second.name;

			fort_printf ("Non-deallocated shared memory found: %s, %d bytes\n",
					name.c_str(), int(size));

			munmap (ptr, size);
		}

		g_shared_mem_map.clear();
	} else {
		//
		// Clean up systemV shared memory

		// Iterate through list of allocations
		list<void*>::iterator iter;
		for (iter = g_shm_list.begin(); iter != g_shm_list.end(); ++iter)
			shmdt(*iter);
		g_shm_list.clear();
	}
#endif
}


// Wrapper to make NAG happy
extern "C" size_t strlen_wrap (const char * str )
{
	return strlen (str);
}



