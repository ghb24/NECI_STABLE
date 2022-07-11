#!/usr/bin/env python
from __future__ import print_function
from functools import reduce
import numpy
import scipy.linalg
import itertools
from pyscf import gto, scf, mcscf
from pyscf.tools import ring
from pyscf.fciqmcscf.fciqmc import *
from pyscf.fci import cistring
from pyscf.fci import rdm
from pyscf.fci import direct_spin1
from pyscf.fci import spin_op
from pyscf.fci import addons

# JPCA, 118, 9925, eq 10

def make_dm12(ci0, norb, nelec):
    def des(ci0, norb, nelec, ap_id, spin):
        if spin == 0:
            ci1 = addons.des_a(ci0, norb, nelec, ap_id)
            ne = nelec[0]-1, nelec[1]
        else:
            ci1 = addons.des_b(ci0, norb, nelec, ap_id)
            ne = nelec[0], nelec[1]-1
        return ci1, ne
    def cre(ci0, norb, nelec, ap_id, spin):
        if spin == 0:
            ci1 = addons.cre_a(ci0, norb, nelec, ap_id)
            ne = nelec[0]+1, nelec[1]
        else:
            ci1 = addons.cre_b(ci0, norb, nelec, ap_id)
            ne = nelec[0], nelec[1]+1
        return ci1, ne

    dm1 = numpy.zeros((norb,2,norb,2))
    for i in range(norb):
        for j in range(norb):
            for i1 in range(2):
                for j1 in range(2):
                    if i1 == j1:
                        ne = nelec
                        ci1, ne = des(ci0, norb, ne, j, j1)
                        ci1, ne = cre(ci1, norb, ne, i, i1)
                        dm1[i,i1,j,j1] = numpy.dot(ci0.ravel(), ci1.ravel())

    dm2 = numpy.zeros((norb,2,norb,2,norb,2,norb,2))
    for i in range(norb):
        for j in range(norb):
            for k in range(norb):
                for l in range(norb):
                    for i1 in range(2):
                        for j1 in range(2):
                            for k1 in range(2):
                                for l1 in range(2):
                                    if i1 + j1 == k1 + l1:
                                        ci1, ne = ci0, nelec
                                        ci1, ne = des(ci1, norb, ne, k, k1)
                                        ci1, ne = des(ci1, norb, ne, l, l1)
                                        ci1, ne = cre(ci1, norb, ne, j, j1)
                                        ci1, ne = cre(ci1, norb, ne, i, i1)
                                        dm2[i,i1,j,j1,k,k1,l,l1] = numpy.dot(ci0.ravel(), ci1.ravel())
    if 0:
        dm1a = numpy.einsum('iajb->ij', dm1)
        dm2a = numpy.einsum('iajbkalb->ijkl', dm2)
        print(abs(numpy.einsum('ipjp->ij', dm2a)/(sum(nelec)-1) - dm1a).sum())
        (dm1a, dm1b), (dm2aa, dm2ab, dm2bb) = \
                direct_spin1.make_rdm12s(ci0, norb, nelec)
        print(abs(dm1a - dm1[:,0,:,0]).sum())
        print(abs(dm2aa - dm2[:,0,:,0,:,0,:,0].transpose(0,2,1,3)).sum())
        print(abs(dm2ab - dm2[:,0,:,1,:,0,:,1].transpose(0,2,1,3)).sum())
        print(abs(dm2ab.transpose(2,3,0,1) - dm2[:,1,:,0,:,1,:,0].transpose(0,2,1,3)).sum())
        print(abs(dm2bb - dm2[:,1,:,1,:,1,:,1].transpose(0,2,1,3)).sum())
        dm2baab = spin_op.make_rdm2_baab(ci0, norb, nelec)
        dm2abba = spin_op.make_rdm2_abba(ci0, norb, nelec)
        print(abs(dm2baab - dm2[:,1,:,0,:,0,:,1].transpose(0,2,1,3)).sum())
        print(abs(dm2abba - dm2[:,0,:,1,:,1,:,0].transpose(0,2,1,3)).sum())
    return dm1, dm2


# evaluate all on orthogonal basis
def dbg_ss_frac(dm1, dm2, norb, mo_coeff, ovlp):
    sinv = numpy.eye(norb)  # on orthogonal basis
    # Note JPCA, 118, 9925, eq 8 may be wrong
    dm2tilde = dm2 + numpy.einsum('jbkc,il,ad->iajbkcld', dm1, sinv, numpy.eye(2))

    sigmax = numpy.zeros((2,2))
    sigmay = numpy.zeros((2,2), dtype=numpy.complex128)
    sigmaz = numpy.zeros((2,2))
    sigmax[0,1] = sigmax[1,0] = 1
    sigmay[0,1] = -1j; sigmay[1,0] = 1j
    sigmaz[0,0] = 1; sigmaz[1,1] = -1
    sigma = numpy.array((sigmax,sigmay,sigmaz))

    sdots = numpy.einsum('xca,xdb->abcd', sigma, sigma).real * .25
    d2 = numpy.einsum('abcd,iajbkcld->ijkl', sdots, dm2tilde)

    def eval_ss_frac(range_A, range_B):
        s_Ax = numpy.zeros_like(ovlp); s_Ax[range_A] = ovlp[range_A]
        s_Bx = numpy.zeros_like(ovlp); s_Bx[range_B] = ovlp[range_B]
        s_Ax = reduce(numpy.dot, (mo_coeff.T, s_Ax, mo_coeff))
        s_Bx = reduce(numpy.dot, (mo_coeff.T, s_Bx, mo_coeff))
        s_xA = s_Ax.T
        s_xB = s_Bx.T

        val = numpy.einsum('ki,lj,ijkl', s_Ax, s_Bx, d2)
        val+= numpy.einsum('ki,lj,ijkl', s_xA, s_xB, d2)
        return val * .5
    return eval_ss_frac

def opt_mag_frac(dm1a, dm1b, norb, nelec, mo_coeff, ovlp):

    sigmax = numpy.zeros((2,2))
    sigmay = numpy.zeros((2,2), dtype=numpy.complex128)
    sigmaz = numpy.zeros((2,2))
    sigmax[0,1] = sigmax[1,0] = 1
    sigmay[0,1] = -1j; sigmay[1,0] = 1j
    sigmaz[0,0] = 1; sigmaz[1,1] = -1
    sigma = numpy.array((sigmax,sigmay,sigmaz))
    # Convert DM into AO basis
    dm_a = reduce(numpy.dot, (mo_coeff, dm1a, mo_coeff.T))
    dm_b = reduce(numpy.dot, (mo_coeff, dm1b, mo_coeff.T))

    def eval_mag_frac(range_A):
        s_Ax = numpy.zeros_like(ovlp)
        # This puts the rows corresponding to the range_A from ovlp into s_Ax
        # First index truncated to be in range_A
        s_Ax[range_A] = ovlp[range_A]
        mag = numpy.zeros((3))
        for spin in range(3):
            mag[spin] = 0.25*(sigma[spin][0,0]*(s_Ax * dm_a.T).sum() + sigma[spin][1,1]*(s_Ax * dm_b.T).sum())
            mag[spin] += 0.25*(sigma[spin][0,0]*(s_Ax.T * dm_a.T).sum() + sigma[spin][1,1]*(s_Ax.T * dm_b.T).sum())
        return mag

    return eval_mag_frac

def opt_ss_frac(dm1a, dm1b, dm2aa, dm2ab, dm2bb, dm2baab, dm2abba, norb, nelec, mo_coeff, ovlp):

    def _bi_trace(dm2, ovlp1, ovlp2):
        return numpy.einsum('jilk,ij,kl->', dm2, ovlp1, ovlp2)

    def eval_ss_frac(range_A, range_B):
        s_Ax = numpy.zeros_like(ovlp); s_Ax[range_A] = ovlp[range_A]
        s_Bx = numpy.zeros_like(ovlp); s_Bx[range_B] = ovlp[range_B]
        s_Ax = reduce(numpy.dot, (mo_coeff.T, s_Ax, mo_coeff))
        s_Bx = reduce(numpy.dot, (mo_coeff.T, s_Bx, mo_coeff))
        s_xA = s_Ax.T
        s_xB = s_Bx.T

        ssz =(_bi_trace(dm2aa, s_Ax, s_Bx)
            - _bi_trace(dm2ab, s_Ax, s_Bx)
            + _bi_trace(dm2bb, s_Ax, s_Bx)
            - _bi_trace(dm2ab.transpose(2,3,0,1), s_Ax, s_Bx)) * .25
        ssz+=(_bi_trace(dm2aa, s_xA, s_xB)
            - _bi_trace(dm2ab, s_xA, s_xB)
            + _bi_trace(dm2bb, s_xA, s_xB)
            - _bi_trace(dm2ab.transpose(2,3,0,1), s_xA, s_xB)) * .25
        ssxy =(_bi_trace(dm2abba, s_Ax, s_Bx)
             + _bi_trace(dm2baab, s_Ax, s_Bx)) * .5
        ssxy+=(_bi_trace(dm2abba, s_xA, s_xB)
             + _bi_trace(dm2baab, s_xA, s_xB)) * .5
        ss = ssxy + ssz
        return ss*.5
    return eval_ss_frac

if __name__ == '__main__':

    mol = gto.Mole()
    mol.atom = [('C', x) for x in ring.make(13,1.25)]
#    mol.atom = '''
#    C  0.00000   0.00000   0.00000
#    C  0.00000   0.00000   1.25000
#    C  0.00000   0.00000   2.5000
#    C  0.00000   0.00000   3.7500
#    '''
    mol.basis='cc-pVDZ'
    mol.symmetry = 1
    mol.symmetry_subgroup = 'C2v'
    mol.charge = 0
    mol.spin = 0
    mol.verbose = 5
    mol.build()
    nbas_atm = 14

    mf = scf.RHF(mol)
    mf.scf()
    mf.analyze()
    
    mc = mcscf.CASCI(mf, 6, 6)

    mo_cas = mf.mo_coeff[:,mc.ncore: mc.ncore + mc.ncas]
#    print 'cas space: ',mo_cas.shape
    norb = mo_cas.shape[1] 
#    print 'norb: ',norb
    ovlp = mc._scf.get_ovlp()
    assert(numpy.allclose(mol.intor('cint1e_ovlp_sph'),ovlp))
    nelec = mc.nelecas

    Addcore = True 
    DoNECI = True 
    
    if DoNECI: 
        mc.fcisolver = FCIQMCCI(mol)
#    mc.canonicalization = False
        mc.fcisolver.generate_neci_input = False

#    emc, e_cas, fcivec = mc.kernel()[:3]
        emc = mc.kernel()[0]
#    ss_tot = spin_op.spin_square(fcivec, norb, nelec)[0]

        # dump orbitals, in case we want to do a calculation having read them in later.
        numpy.save('neci-orbs',mf.mo_coeff)

#    ss_tot = spin_op.spin_square(fcivec, norb, nelec)[0]
        ss_tot = 0.0

        reorder=False
        dm1a, dm1b = read_neci_1dms(mc.fcisolver,norb,nelec)
        dm2aa, dm2ab, dm2bb = read_neci_2dms(mc.fcisolver,norb,nelec,reorder=reorder)
        if Addcore:
            # Add core
            dm1a, dm1b, dm2aa, dm2ab, dm2bb = add_spinned_core_rdms(mf, mc.ncore, dm1a, dm1b, dm2aa, dm2ab, dm2bb,reorder=reorder)
            nelec_full = (nelec[0]+mc.ncore, nelec[1]+mc.ncore)
            occorbs = mf.mo_coeff[:,:mc.ncore + norb]
            n_occorbs = occorbs.shape[1]
        else:
            nelec_full = nelec
            occorbs = mo_cas
            n_occorbs = occorbs.shape[1]

        dm2abba = dm2ab.copy()
        dm2abba = -dm2abba.transpose(2,1,0,3)
        # This gives the 'normal' baab matrix. However, this calls the underscore one (i.e. no reordering)
        for i in range(dm2abba.shape[0]):
            dm2abba[:,i,i,:] += dm1b 
        dm2baab = dm2abba

        if reorder:
            spinfree_2 = dm2aa + dm2bb + 2*dm2ab
            one_sf, two_sf = find_full_casscf_12rdm(mc.fcisolver, mf.mo_coeff, 'spinfree_TwoRDM.1', norb, nelec)
            assert(numpy.allclose(two_sf[:n_occorbs,:n_occorbs,:n_occorbs,:n_occorbs], spinfree_2))
    
    mc2 = mcscf.CASCI(mf, 6, 6)
#    mc2.canonicalization = False
    mc2.fcisolver = direct_spin1.FCISolver(mol)
    emc2, e_cas2, fcivec2 = mc2.kernel()[:3]
    #ss_tot = spin_op.spin_square(fcivec2, norb, nelec)[0]
#    print 'ss_tot exact: ',ss_tot
    (dm1a_, dm1b_), (dm2aa_, dm2ab_, dm2bb_) = direct_spin1.make_rdm12s(fcivec2, norb, nelec, reorder=False)
    
    #f_dbg = dbg_ss_frac(dm1, dm2, norb, mo_cas, ovlp)
    # dm2baab(pq,rs) = <p(beta)* q(alpha) r(alpha)* s(beta)
    dm2baab_ = spin_op._make_rdm2_baab(fcivec2, norb, nelec)
    # dm2abba(pq,rs) = <q(alpha)* p(beta) s(beta)* r(alpha)
    dm2abba_ = spin_op._make_rdm2_abba(fcivec2, norb, nelec)
    ss_tot = spin_op.spin_square(fcivec2, norb, nelec)[0]

    f_opt_act = opt_ss_frac(dm1a_, dm1b_, dm2aa_, dm2ab_, dm2bb_, dm2baab_,  
            dm2abba_, norb, nelec, mo_cas, ovlp)
    f_mag_act = opt_mag_frac(dm1a_, dm1b_, norb, nelec, mo_cas, ovlp)
    if DoNECI:
        f_opt_full = opt_ss_frac(dm1a, dm1b, dm2aa, dm2ab, dm2bb, dm2baab,  
                dm2abba, n_occorbs, nelec_full, occorbs, ovlp)
        f_mag_full = opt_mag_frac(dm1a, dm1b, n_occorbs, nelec_full, occorbs, ovlp)
    if not Addcore and DoNECI:
        assert(numpy.allclose(dm1a,dm1a_))
        assert(numpy.allclose(dm1b,dm1b_))
        assert(numpy.allclose(dm2aa,dm2aa_))
        assert(numpy.allclose(dm2ab,dm2ab_))
        assert(numpy.allclose(dm2bb,dm2bb_))
        assert(numpy.allclose(dm2baab,dm2baab_))
        assert(numpy.allclose(dm2abba,dm2abba_))
        assert(numpy.allclose(mo_cas,occorbs))

    # natm atoms, nbas_atm AO for each
    ss_opt_full = 0.0
    ss_opt_act = 0.0
    mag_opt_act = [0.0, 0.0, 0.0]
    mag_opt_full = [0.0, 0.0, 0.0]
    magz_tot = 0.0
    iatm = 0
    for iA in range(0, mo_cas.shape[0], nbas_atm):
        mag_opt_act += f_mag_act(range(iA,iA+nbas_atm))
        print('active space magnetization on atm {} is {} '.format(iatm,f_mag_act(range(iA,iA+nbas_atm))))
        if DoNECI: mag_opt_full += f_mag_full(range(iA,iA+nbas_atm))
        iatm += 1
        for iB in range(0, mo_cas.shape[0], nbas_atm):
#            ss_dbg += f_dbg(range(iA,iA+nbas_atm), range(iB,iB+nbas_atm))
            ss_opt_act += f_opt_act(range(iA,iA+nbas_atm), range(iB,iB+nbas_atm))
            if DoNECI: ss_opt_full += f_opt_full(range(iA,iA+nbas_atm), range(iB,iB+nbas_atm))
    print('Total S^2 of the system (tot, act, full): ',ss_tot, ss_opt_act, ss_opt_full)
    print('Total Sz of the system (tot, act, full): ',magz_tot, mag_opt_act[0], mag_opt_full[0])

    #Pick first atom, and examine the spin correlation functions to neighbouring atoms
    iatm = 0
    spincorr_act = []
    spincorr_full = []
    for iB in range(0, mo_cas.shape[0], nbas_atm):
        if iatm == 0: 
            probe_act = f_opt_act(range(nbas_atm), range(iB,iB+nbas_atm))
            if DoNECI: probe_full = f_opt_full(range(nbas_atm), range(iB,iB+nbas_atm))
        spincorr_act.append(f_opt_act(range(nbas_atm), range(iB,iB+nbas_atm)))
        print('active space spin-corr func to atm {} is {}   {}'.format(iatm,spincorr_act[-1],spincorr_act[-1]/spincorr_act[0]))
        if DoNECI: 
            spincorr_full.append(f_opt_full(range(nbas_atm), range(iB,iB+nbas_atm)))
            print('full space spin-corr func to atm {} is {}   {}'.format(iatm,spincorr_full[-1],spincorr_full[-1]/spincorr_full[0]))
        iatm += 1

    print('')
    print()
    if len(spincorr_act) == iatm:
        print('Active space spin-correlation function: ')
        print('Atom         SS-function       Normalized SS-func')
        for i in range(iatm):
            print('{}            {}             {}'.format(i,spincorr_act[i],spincorr_act[i]/spincorr_act[0]))
    print('')
    print()
    if len(spincorr_full) == iatm:
        print('Full space spin-correlation function: ')
        print('Atom         SS-function       Normalized SS-func')
        for i in range(iatm):
            print('{}            {}             {}'.format(i,spincorr_full[i],spincorr_full[i]/spincorr_full[0]))

    with open('spincorrfns','w') as f:
        f.write('# Atom  SS-act  SS-act-norm  SS-full  SS-full-norm\n')
        for i in range(iatm):
            f.write('{}     {}     {}    {}    {}\n'.format(i,spincorr_act[i],spincorr_act[i]/spincorr_act[0],spincorr_full[i],spincorr_full[i]/spincorr_full[0]))
