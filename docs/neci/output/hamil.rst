.. _output_hamil:

-----
HAMIL
-----

HAMIL contains the non-zero elements of the Hamiltonian in three columns::

    i  j  Hij

where: 

    Hij :math:`=\bra D_{\veci} | H | D_{\vecj} \ket`.

    i,j are indices for the :math:`\veci,\vecj` determinants, with i:math:`\le`
    j, and increasing i.

    :math:`D_i` corresponds to the i-th determinant as given in DETS.


