#include "lib/types.h"

// Fortran requires stdc
extern "C" {

// The fortran main function
void performfcimcycpar (void (*)());

// Annihilation
void (*g_annihilate)() = 0;

void set_annihilator (void (*annihilate)())
{
	g_annihilate = annihilate;
}

void call_fcimc_cyc_par ()
{
	if (!g_annihilate)
		stop_all (__FUNCTION__, "Annihilation routine not set.");

	performfcimcycpar (g_annihilate);
}

}
