title: Trial Wavefunctions 
---

## Trial wave functions

By default, NECI uses a single reference determinant, \(|D_0\rangle\),
in the projected energy estimator, or potentially a linear combination
of two determinants if the the HPHF code is being used.
\[E_0 = \frac{\langle D_0 | \hat{H} | \Psi\rangle}{\langle D_0 | \Psi}\rangle .\]
This estimator can be improved by using a more accurate estimate of the
true ground state, a trial wave function, \(|\Psi^T\rangle\),
\[E_0 = \frac{\langle{\Psi^T | \hat{H} | \Psi}\rangle}{\langle{\Psi^T | \Psi}\rangle}.\]
Such a trial wave function can be used in NECI using by adding the
`trial-wavefunction` option to the Calc block. You must also specify a
trial space. The trial wave function used will be the ground state of
the Hamiltonian projected into this trial space.

The trial spaces available are the same as the core spaces available for
the semi-stochastic option. However, you must replace `core` with
`trial`. For example, to use all single and double excitations of the
reference determinant, one should use the ‘doubles-trial’ option.
