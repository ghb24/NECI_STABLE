.. _input_examples:

-----------------------
INITIATOR CALCULATIONS
-----------------------

Molecular systems
-----------------

For molecular systems, a valid input file is of the form::

    TITLE

    SYSTEM READ NOORDER
    SYMIGNOREENERGIES
    ELECTRONS 9
    NONUNIFORMRANDEXCITS
    SPIN-RESTRICT -1
    NOBRILLOUINTHEOREM
    ENDSYS

    CALC
    METHODS
    METHOD VERTEX FCIMC
    ENDMETHODS
    BETA 500
    NMCYC 5000000
    DIAGSHIFT 0.D+0
    INITWALKERS 2000000
    TAU 0.00005
    MAXNOATHF 50000 10000
    STARTSINGLEPART
    TRUNCATECAS 5 7
    TRUNCINITIATOR
    KEEPDOUBSPAWNS
    ADDTOINITIATOR 3
    DIRECTANNIHILATION
    ANNIHILATERANGE
    SHIFTDAMP 0.10
    GROWMAXFACTOR 2.0
    CULLFACTOR 2.0
    STEPSSHIFT 20
    GLOBALSHIFT
    MEMORYFACPART 2.0
    MEMORYFACSPAWN 0.5
    REGENEXCITGENS
    RHOEPSILON 0.D+00
    PATHS 1
    ENDCALC

    INTEGRAL
    FREEZE 0,0
    ENDINT

    LOGGING
    BINARYPOPS
    ENDLOG
    END

This is a flourine atom calculation using an initial fixed initiator space of (5,6) (meaning 5 electrons within
6 spatial orbitals).  This is a restricted open shell calculation with 1 unpaired electron (**SPIN-RESTRICT**),
and so brillouins theorem does not hold (**NOBRILLOUINTHEOREM**), although without this keyword the code should 
figure this out on it's own.
This calculation uses the initiator method which is begun by using **TRUNCINITIATOR**.
Determinants are added to the initiator space when their absolute population is greater than 3
(**ADDTOINITIATOR**).
An additional rule is included whereby if two different determinants spawn on the same determinant with the same sign, 
the spawned walkers are allowed to live (**KEEPDOUBSPAWNS**).
This calculation starts with a single walker on the HF determinant (**STARTSINGLEPART**) and allows this population 
to rise to 50,000 and then lets the shift change so that the total number of walkers stays constant at this point.  
However, if the HF population drops more than 10,000 below this fixed value, the shift is held constant again and 
the number of walkers allowed to grow until the HF population again reaches 50,000 (**MAXNOATHF**).
At the end of the calculation a POPSFILE will be printed by default, and in this case this will be in binary 
(**BINARYPOPS**).

