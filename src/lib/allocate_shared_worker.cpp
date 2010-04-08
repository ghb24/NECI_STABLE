#include <stdint.h>
#include <stdarg.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <map>
#include <string>
using std::map;
using std::string;

// Prototype of stop_all function.
extern "C" void stop_all_ (const char* a, const char* b);

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
	map_det_t (string name_, int size_) : name(name_), size(size_) {}
	string name;
	size_t size;
};
map<void*,map_det_t> g_shared_mem_map;

//
// This function acquires a named region of shared memory, of the specified
// size, allocating or resizing it if necessary. It then returns this as a
// standard C pointer --> requires fortran 2003 features to use.
extern "C" void alloc_shared_worker_ (const char * name, void ** ptr,
		                              const size_t size)
{
	// Acquire a named shared memory file descriptor. nb. Enforce 10 character
	// name limit.
	int fd = shm_open (name, O_CREAT | O_RDWR, (mode_t)0600);
	if (fd == -1)
		stop_all_ (__FUNCTION__, "Creating shared memory object failed.");

	// Set the length of the named region to the required length
	if (ftruncate(fd, size) == -1)
		stop_all_ (__FUNCTION__, "Setting size of shared memory regoin failed.");

	// Map the region into the current address space.
	*ptr = mmap (NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, 
			fd, 0);
	if (*ptr == MAP_FAILED)
		stop_all_ (__FUNCTION__, "Mapping shared memory failed.");

	// Once we have mapped the region, it will remain available even when the
	// file descriptor has closed.
	close(fd);

	g_shared_mem_map[*ptr] = map_det_t(name, size);
}

//
// Given a C pointer to a region of shared memory, retrieve all of the
// relevant details from the map, and then deallocate the memory and remove
// the mapping.
extern "C" void dealloc_shared_worker_ (void * ptr)
{
	// Find the shared memory in the global list.
	map<void*,map_det_t>::iterator det;
	det = g_shared_mem_map.find(ptr);

	if (det == g_shared_mem_map.end())
		stop_all_ (__FUNCTION__, "The specified shared memory was not found.");

	munmap (ptr, det->second.size);

	shm_unlink (det->second.name.c_str());

	g_shared_mem_map.erase(det);
}

//
// Clean up any shared allocations which have not been properly deallocated.
extern "C" void cleanup_shared_alloc_ ()
{
	// Iterate through the list of shared allocations
	map<void*,map_det_t>::iterator iter;
	for (iter = g_shared_mem_map.begin(); iter != g_shared_mem_map.end();
	     ++iter) {
		size_t size = iter->second.size;
		void * ptr = iter->first;
		string name = iter->second.name;

		printf ("Non-deallocated shared memory found: %s, %d bytes\n",
				name.c_str(), int(size));

		munmap (ptr, size);
		shm_unlink (name.c_str());
	}

	g_shared_mem_map.clear();
}

