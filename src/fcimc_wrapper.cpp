#include "lib/types.h"

// Fortran requires stdc
extern "C" {

// The fortran main function
void performfcimcycpar (void (*)(), void (*)(), int32_t, void (*)(),
		                void (*)(), void (*)());

//
// A useful function which does nothing, and is therefore a useful default
// for a subroutine, if we sometimes don't want anything.
// cdecl is quite useful, as the caller cleans up the stack. This means we
// don't need to worry about arguments etc.
void null_function () {}

//
// Store the function pointers
static void (*g_annihilate)() = 0;
static void (*g_generate_excitation)() = 0;
static void (*g_get_spawn_helement)() = 0;
static void (*g_attempt_create)() = 0;
static void (*g_encode_child)() = &null_function;
static int32_t g_max_excit_level = 2;

//
// Setter functions for function pointered bits
void set_annihilator (void (*annihilate)())
	{ g_annihilate = annihilate; }
void set_excit_generator (void (*generate_excitation)()) 
	{ g_generate_excitation = generate_excitation; }
void set_get_spawn_helement (void (*get_spawn_helement)())
	{ g_get_spawn_helement = get_spawn_helement; }
void set_max_excit_level (int32_t level)
	{ g_max_excit_level = level; }
void set_attempt_create (void (*attempt_create)())
	{ g_attempt_create = attempt_create; }
void set_encode_child (void (*encode_child)())
	{ g_encode_child = encode_child; }



void call_fcimc_cyc_par ()
{
	if (!g_annihilate)
		stop_all (__FUNCTION__, "No annihilation routine set.");

	if (!g_generate_excitation)
		stop_all (__FUNCTION__, "No excitation generation routine set.");

	if (!g_attempt_create)
		stop_all (__FUNCTION__, "No attempt_create routine set.");

	if (!g_get_spawn_helement)
		stop_all (__FUNCTION__, "No routine to calculate matrix elements for "
				                "excitations is set.");

	performfcimcycpar (g_annihilate, g_generate_excitation, 
	                   g_max_excit_level, g_attempt_create, 
					   g_get_spawn_helement, g_encode_child);
}

}
