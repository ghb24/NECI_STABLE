#include "lib/types.h"

// Fortran requires stdc
extern "C" {

// The fortran main function
void performfcimcycpar (void (*)(), void (*)());

//
// Store the function pointers
static void (*g_annihilate)() = 0;
static void (*g_generate_excitation)() = 0;

//
// Setter function for annihilation
void set_annihilator (void (*annihilate)())
{
	g_annihilate = annihilate;
}

//
// Setter function for excitation generation
void set_excit_generator (void (*generate_excitation)())
{
	g_generate_excitation = generate_excitation;
}


void call_fcimc_cyc_par ()
{
	if (!g_annihilate)
		stop_all (__FUNCTION__, "No annihilation routine set.");

	if (!g_generate_excitation)
		stop_all (__FUNCTION__, "No excitation generation routine set.");

	performfcimcycpar (g_annihilate, g_generate_excitation);
}

}
