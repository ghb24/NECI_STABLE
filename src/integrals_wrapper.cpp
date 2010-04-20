#include "lib/types.h"


// Type of a getumatel call.
typedef helement_t (*fn_getumatel_t)(int32_t i, int32_t j, int32_t k,
		                             int32_t l);

//
// Store function pointer (stack) for getumatel.
fn_getumatel_t g_get_umat_el[2] = {0, 0};

// The array G1 contains basis function information. Use a reference to make
// it easily accessible.
extern basisfn_t * systemdata_mp_g1_;
basisfn_t * &g1 = systemdata_mp_g1_;

//
// __cdecl linkage to interoperate with Fortran.
extern "C" {

void set_getumatel_fn (fn_getumatel_t fn)
{
	// This allows you to create a nifty stack of functions (e.g. for tFixLz)
	if (g_get_umat_el[0])
		g_get_umat_el[1] = g_get_umat_el[0];

	g_get_umat_el[0] = fn;
}

//
// If we store spin orbitals and are fixing Lz, then do a quick check to see
// if we know that umatel = 0. Then call the next function on the stack.
helement_t get_umat_el_fixlz_storespinorbs (int32_t i, int32_t j, int32_t k,
                                            int32_t l)
{
	// If ew are fixing Lz, then <ij|kl> != <kj|il> necessarily, since we
	// have complex orbitals (though real integrals) and want to ensure that
	// we conserve momentum. i.e. momentum of bra = mom of ket.
	if ( (g1[i-1].Ml + g1[j-1].Ml) != (g1[k-1].Ml + g1[l-1].Ml) )
		return 0;

	return g_get_umat_el[1](i, j, k, l);
}

//
// If we are not storing spin orbitals, and are fixing Lz, then as above
helement_t get_umat_el_fixlz_notspinorbs (int32_t i, int32_t j, int32_t k,
                                          int32_t l)
{
	if ( (g1[2*i-1].Ml + g1[2*j-1].Ml) != (g1[2*k-1].Ml + g1[2*l-1].Ml) )
		return 0;

	return g_get_umat_el[1](i, j, k, l);
}

helement_t get_umat_el (int32_t i, int32_t j, int32_t k, int32_t l)
{
	return g_get_umat_el[0](i, j, k, l);
}


}
