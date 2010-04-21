#include "lib/types.h"

// Fortran requires stdc
extern "C" {

// The fortran main function
void performfcimcycpar (void (*)(), void (*)(), void (*)(), void (*)(),
                        void (*)());

//
// A useful function which does nothing, and is therefore a useful default
// for a subroutine, if we sometimes don't want anything.
// cdecl is quite useful, as the caller cleans up the stack. This means we
// don't need to worry about arguments etc.
void null_function () {}

//
// Store the function pointers
//
// ***************** NOTE and WARNING ********************
// These function pointers must NEVER be called directly from C. They are
// intentionally declared with incorrect argument lists and return values 
// (i.e. none), as they are ONLY to be used from Fortran.
//
// --> The correct parameter etc. can be found in the interfaces for the
//     set_* functions declared at the top of FciMCPar.F90.
// --> They must only be used by being passed as arguments to a fortran
//     routine with _matching_ interfaces --> can then be called from there
//
// This is a slight abuse of __cdecl, which does not enforce type details
// when linking --> only the correct call address, and function name need to
// be present in files for linking. Therfore we can cheat across language
// boundaries.
// *******************************************************

static void (*g_generate_excitation)() = 0;
static void (*g_get_spawn_helement)() = 0;
static void (*g_attempt_create)() = 0;
static void (*g_encode_child)() = &null_function;
static void (*g_new_child_stats)() = 0;

//
// Setter functions for function pointered bits
void set_excit_generator (void (*generate_excitation)()) 
	{ g_generate_excitation = generate_excitation; }
void set_get_spawn_helement (void (*get_spawn_helement)())
	{ g_get_spawn_helement = get_spawn_helement; }
void set_attempt_create (void (*attempt_create)())
	{ g_attempt_create = attempt_create; }
void set_encode_child (void (*encode_child)())
	{ g_encode_child = encode_child; }
void set_new_child_stats (void (*new_child_stats)())
	{ g_new_child_stats = new_child_stats; }



void call_fcimc_cyc_par ()
{
	if (!g_generate_excitation)
		stop_all (__FUNCTION__, "No excitation generation routine set.");

	if (!g_attempt_create)
		stop_all (__FUNCTION__, "No attempt_create routine set.");

	if (!g_get_spawn_helement)
		stop_all (__FUNCTION__, "No routine to calculate matrix elements for "
				                "excitations is set.");

	if (!g_new_child_stats)
		stop_all (__FUNCTION__, "No routine for new child statistics is set.");
#ifdef PARALLEL
	performfcimcycpar (g_generate_excitation, g_attempt_create, 
	                   g_get_spawn_helement, g_encode_child,
					   g_new_child_stats);
#endif
}

}
