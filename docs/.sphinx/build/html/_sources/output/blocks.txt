.. _output_blocks:

------
BLOCKS
------

BLOCKS contains the list of blocks used if the **BLOCK** option is used to
calculate the Hamiltonian in the form::

    BlockIndex   K(1)   K(2)   K(3)   MS SYM   nDets nTotDets

where:

    BlockIndex starts from 1.

    K(1:3) are momentum values for the **HUBBARD** and **UEG** symmetries and
    are irrelevant for other systems.

    MS is :math:`2S_z`.

    Sym is the spatial symmetry index of the determinant.  It can be a single
    number or a set of numbers enclosed in parentheses.  Relevant for systems
    other than **HUBBARD** and **UEG**.

    nDets is the number of symmetry unique determinants in the block.

    nTotDets is the total number of determinants in the block.

