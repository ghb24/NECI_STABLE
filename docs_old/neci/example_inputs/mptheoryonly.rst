.. _input_examples:

---------------------------
Standalone MP2 calculations
---------------------------

The **CALC** block is simply::
   
    Calc
        MPTheory only
    EndCalc

CPMD
----

For a straight-forward **CPMD**-based calculation using, say, 200 spin-virtuals, the input file is::

    System CPMD
    EndSys

    Calc
        MPTheory only
    EndCalc

    Integrals
        UMatCache MB 750
        Freeze 0,-200
    EndInt

where we have also allocatabed a maximum of 750MB to be used in caching the integrals.

If we wish to ignore the single excitations of the reference (Kohn--Sham) determinant, then we can use::

    System CPMD
    EndSys

    Calc
        MPTheory only
        Excitations doubles
    EndCalc

    Integrals
        UMatCache MB 750
        Freeze 0,-200
    EndInt

Molecular systems
-----------------

Similarly, for molecular systems, a valid input file is of the form::

    System READ
        Electrons 36
    EndSys

    Calc
        MPTheory only
    EndCalc
