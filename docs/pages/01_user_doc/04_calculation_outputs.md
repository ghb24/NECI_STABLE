---
title: Output files
---


## Output files

Apart from the output that is printed to standard out
  there are some other useful files created by `NECI`.
Some of them are only created, when certain keywords are given.

### FCIMCStats

This file contains whitespace delimited data that is written
  every 10nth iteration.
Currently (2021-06-08) there are 35 columns.
The information in this file is useful in virtually every way of using `NECI`.

The columns are listed in the following.
Important and often used columns are highlighted red,
rarely used columns are greyed out.

**\textcolor{red}{1 Steps (Step)}**<br>
The number of iterations.

**\textcolor{red}{2 Shift (Shift)}**<br>
The shift \(S\) for population control. Equals the correlation
energy in the equilibrium.
Shift + Reference energy should equal the projected energy in equilibrium.

**3 Walker Change (WalkerCng)**<br>
Absolute change of the walker number.
(Since there are fractional walkers, this can be non-integer.)

**4 Growth Rate (GrowRate)**<br>
Relative change of walker number.

**\textcolor{red}{5 Total Walkers (TotWalkers)}**<br>
The number of total walkers.

**6 Annihilation (Annihil)**<br>
The number of annihilated walkers.

**\textcolor{gray}{7 Number of died walkers (NoDied)}**<br>
\textcolor{gray}{The number of died walkers (from annihilation or diagonal death step).}

**\textcolor{gray}{8 Number of born walkers (NoBorn)}**<br>
\textcolor{gray}{The number of born walkers.}

**\textcolor{gray}{9 Projected averaged correlation energy (Proj.E)}**<br>
\textcolor{gray}{Averaged correlation energy as calculated from the projected energy expression.  The unaveraged values are given in column 11.}

**\textcolor{gray}{10 Average Shift (Av.Shift)}**<br>
\textcolor{gray}{Averaged shift.  Note that it should not be used, when the average was taken over a period that not in stationary conditions.}


**\textcolor{red}{11 Projected instantaneous correlation energy (Proj.E.ThisCyc)}**<br>
Instantaneous (not averaged) correlation energy as calculated
from the projected energy expression.
\begin{equation}
    E_{\text{proj}} - H_{00}
  =
    \frac{\left\langle D_0 | H | \Psi \right\rangle} {\left\langle D_0  \Psi \right\rangle} - H_{00}
  =
    \sum_{j\neq 0} H_{0j} \frac{C_j}{C_0} \quad.
\end{equation}
The averaged values are given in column 9.
The total energy (correlation + reference energy) is given in column 23)

**\textcolor{red}{12 Number of walkers at reference determinant \(D_0\) (NoatHF)}**<br>
The number of walkers at the reference determinant.
(Since there are fractional walkers, this can be non-integer.)

<!-- Unimportant -->
**\textcolor{gray}{13 Number of walkers at doubles (NoatDoubs)}**<br>
\textcolor{gray}{Number of walkers that are a double excitation away from the reference.}

**14 Acceptance rate (AccRat)**<br>
Probability that a spawn gets accepted.
This is **not** the condition probability, but the absolute one.

**15 Uniqe determinants (UniqueDets)**<br>
Number of unique configurations (determinants/CSFs).

**16 Iteration time (IterTime)**<br>
This is the time averaged over the last 10 iterations.

<!-- Unimportant -->
**\textcolor{gray}{17 (FracSpawnFromSing)}**<br>
TODO

<!-- Unimportant -->
**\textcolor{gray}{18 (WalkersDiffProc)}**<br>
TODO

**\textcolor{red}{19 Total imaginary time \(\tau\) (TotImagTime)}**<br>
The elapsed **imaginary** time since start of the dynamics.

<!-- Confirm -->
<!-- Broken -->
**20 (ProjE.ThisIter)**<br>
TODO

<!-- Unimportant -->
**21 (HFInstShift)**<br>
TODO
<!-- Ask Khaldoon -->

<!-- Unimportant -->
**22 (TotInstShift)**<br>
TODO
<!-- Ask Khaldoon -->

**\textcolor{red}{23 Projected instantaneous total energy (Tot-Proj.E.ThisCyc)}**<br>
Instantaneous (not averaged) total energy as calculated
from the projected energy expression
\begin{equation}
    E_{\text{proj}}
  =
    \frac{\left\langle D_0 | H | \Psi \right\rangle} {\left\langle D_0  \Psi \right\rangle}
  =
    H_{00} + \sum_{j\neq 0} H_{0j} \frac{C_j}{C_0} \quad.
\end{equation}
The instantaneous correlation energy is given in column 11.
The nominator and denominator are given in columns 24 and 25 and should be
used for statistical analysis of errors.
It is averaged over the last 10 iterations.

**24 Instantaneous denominator of projected energy (HFContribtoE)**<br>
This is the instantaneous denominator of the projected energy \(C_0\)
  and equivalent to the reference weight (column 12).
Column 24 and 25 are used to evaluate the other projected energy columns
  (9, 11, and 23).
It is averaged over the last 10 iterations.

**25 Instantaneous numerator of projected energy (NumContribtoE)**<br>
This is the instantaneous nominator of the projected energy
  \begin{equation}
    \sum_{j\neq 0} H_{0j} C_j
  \end{equation}
Column 24 and 25 are used to evaluate the other projected energy columns
  (9, 11, and 23).
It is averaged over the last 10 iterations.


<!-- Confirm -->
**\textcolor{red}{26 Amplitude of reference (HF weight)}**<br>
Unlike the name suggests, this is the **amplitude** weight of the reference determinant
  given by:
  \begin{equation}
      \frac{|c_0|}{|\Psi|}
    =
      \frac{|c_0|}{\sqrt{\sum_i |c_i|^2}}
    \quad.
  \end{equation}
This means that this column can be obtained by dividing column 12 by column 27.

**27 \(| \Psi |\) (|Psi|) **<br>
Instantaneous L2-norm of the current wavefunction.
\begin{equation}
  | \Psi | = \sqrt{\sum_i |c_i|^2}
\end{equation}

**28 Expectation value of \(S^2\) operator over all determinants (Inst S^2)**<br>
Requires the `instant-s2-full` keyword in the logging block.

**29 Expectation value of \(S^2\) operator over initiator determinants (Inst S^2)**<br>
Requires the `instant-s2-init` keyword in the logging block.

<!-- Unimportant -->
<!-- Confirm -->
**\textcolor{gray}{30 Absolute instantaneous projected correlation energy (AbsProjE)}**<br>
\textcolor{gray}{L1-norm of projected correlation energy expression.}
\begin{equation}
  \sum_{j\neq 0} |H_{0j}| \frac{|C_j|}{|C_0|} \quad.
\end{equation}

<!-- Unimportant -->
<!-- Uninitialised Garbage -->
**\textcolor{gray}{31 (PartsDiffProc)}**<br>
\textcolor{gray}{Uninitialised Garbage.}

**32 Weight of semistochastic space (|Semistoch|/|Psi|)**<br>
If \(S_C\) is the semistochastic or core space, then we have:
\begin{equation}
  \frac{\sum_{j \in S_C} |c_j|^2} {\sum_i |c_i|^2} \quad.
\end{equation}

<!-- Unimportant -->
**\textcolor{gray}{33 Largest spawn per iteration (MaxCycSpawn)}**<br>
\textcolor{gray}{This is the largest spawn per iteration.}

**34 The number of discarded excitations (InvalidExcits)**<br>
This is the number of discarded excitations averaged over the last 10 iterations.

**35 The number of valid excitations (ValidExcits)**<br>
This is the number of valid excitations averaged over the last 10 iterations.
