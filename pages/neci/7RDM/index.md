title: RDM Generation 
---

## RDM generation

Currently the 2-RDMs can only be calculated for closed shell systems.
However, calculation and diagonalisation of only the 1-RDM is set up for
either open shell or closed shell systems.

The original theory behind the calculation of the RDMs (including
details of parallelisation) can be found in the paper:
<http://arxiv.org/abs/1410.6047>. The most accurate RDM method (which is
also unbiased) is the double-run approach, which requires the code to be
compiled with the `-D__DOUBLERUN` flag in `CPPFLAGS` in the Makefile.
This propagates two completely independent populations of walkers, and
calculates an unbiased RDM by taking cross-terms between the two
populations.

The calculation of the diagonal elements is done by keeping track of the
average walker populations of each occupied determinant, and how long it
has been occupied. The diagonal element from Di is then calculated as
&lt;Ni&gt;(pop1) x &lt;Ni&gt;(pop2) x \[No. of iterations in this
block\], and this is included every time we start a new averaging block,
which can occur when a determinant becomes unoccupied in either
population, or when we require the calculation of the RDM energy during
the simulation. As such, the exact RDM accumulated is dependent on the
interval of RDM energy calculations, but in an unbiased way.

The off diagonal elements are sampled through spawning events, and use
instantaneous walker populations. Subtle details in the code are:

1\. RDMs take contributions to diagonal elements and HF connections
whenever the RDM energy is calculated. As the averaging blocks are now
reset at this point too, this change is unbiased. 2. The off-diagonal
contributions (both HF connections and other contributions sampled
through spawning events) contain contributions from both cross-terms.
I.e \(N_i(pop1)*N_j(pop2)\) as well as \(N_i(pop2)*N_j(pop1)\).

There is also the facility to do a single-run calculation of the RDM.
This method is BIASED, so should not be used for high accuracy
calculations. However, it is cheaper than the double run method, both in
memory and in simulation time, so may be useful for rough-and-ready
calculations. Reducing the effect of the bias in the SR method can be
done by applying a cutoff to the diagonal contributions, such that
contributions are only added in if the average sign of the determinant
at the time of adding in the contribution exceeds some preset parameter.
When calculating RDMs, the RDM energy will be printed at the end of the
calculation, which is one measure of the accuracy of the RDMs. Also
printed by default are the maximum error in the hermiticity
(2-RDM(i,j;a,b) - 2-RDM(a,b;i,j)) and the sum of the absolute errors.

### Reading in / Writing out the RDMs for restarting calculations

Two types of 2-RDMs can be printed out. The final normalised hermitian
2-RDMs of the form `TwoRDM_a***`, or the binary files
`TwoRDM_POPS_a***`, which are the unnormalised RDMs, before hermiticity
has been enforced. The first are the `2-RDM(i,j;a,b)` matrices, which
are printed in spatial orbitals with \(i<j\), \(a<b\) and \(i,j<a,b\).
The second are the ones to read back in if a calculation is restarted
(they are also printed in spatial orbitals with \(i<j\) and \(a<b\), but
for both \(i,j,a,b\) and \(a,b,i,j\) because they are not yet
hermitian). These are the matrices exactly as they are at that point in
the calculation. By default the final normalised 2-RDMs will always be
printed, and the `TwoRDM_POPS_a***` files are connected to the
`popsfile/binarypops` keywords - i.e. if a wavefunction popsfile is
being printed and the RDMs are being filled, a RDM `POPSFILE` will be
also. If only the 1-RDM is being calculated, `OneRDM_POPS/OneRDM` files
will be printed in the same way.
