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

The columns are:

<!-- # 28.Inst S^2 29.Inst S^2   30.AbsProjE   31.PartsDiffProc 32.|Semistoch|/|Psi|  33.MaxCycSpawn   34.InvalidExcits  35. ValidExcits -->
**1 Steps**<br>
The number of iterations.

**2 Shift**<br>
The shift \(S\) for population control. Equals the correlation
energy in the equilibrium.
Shift + Reference energy should equal the projected energy in equilibrium.

**3 Walker Change**<br>
TODO

**4 Growth Rate**<br>
TODO

**5 Total Walkers**<br>
The number of total walkers.

**6 Annihilation**<br>
The number of annihilated walkers.

**7 Number of died walkers**<br>
The number of died walkers (from annihilation **or** diagonal death step).

**8 Number of born walkers**<br>
The number of born walkers.

**9 Projected energy**<br>
TODO

**10 Average Shift**<br>
TODO

**11 Projected energy of this cycle**<br>
TODO

**12 Number of walkers at reference determinant \(D_0\)**<br>
The number of walkers at the reference determinant.

**13 NoatDoubs**<br>
TODO

**14 AccRat**<br>
TODO

**15 UniqueDets**<br>
TODO

**16 Iteration time**<br>
This is the time averaged over the last 10 iterations.

**17 FracSpawnFromSing**<br>
TODO

**18 WalkersDiffProc**<br>
TODO

**19 Total imaginary time**<br>
TODO

**20 ProjE.ThisIter**<br>
TODO

**21 HFInstShift**<br>
TODO

**22 TotInstShift**<br>
TODO

**23 Tot-Proj.E.ThisCyc**<br>
TODO

**24 HFContribtoE**<br>
TODO

**25 NumContribtoE**<br>
TODO

**26 HF weight**<br>
TODO

**27 \(| \Psi |\)**<br>
TODO

**28 Inst S^2**<br>
TODO

**29 Inst S^2**<br>
TODO

**30 AbsProjE**<br>
TODO

**31 PartsDiffProc**<br>
TODO

**32 |Semistoch|/|Psi|**<br>
TODO

**33 MaxCycSpawn**<br>
TODO

**34 InvalidExcits**<br>
TODO

**35 ValidExcits**<br>
TODO


