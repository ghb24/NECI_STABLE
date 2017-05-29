.. _output_classpaths2:

-----------
CLASSPATHS2
-----------

CLASSPATHS2 is calculated when a vertex sum is performed and gives a histogram of the weights for each graph type.

For each graph type, a header is printed out followed by the histogram data.

Header::

    Class nGs   TotWeightPos   TotWeightNeg

Body::

    log_10(weight)   nPos   nNeg

where:

    Class, nGs, TotWeightPos, TotWeightNeg are the same as in CLASSPATHS.

    The body consists of lines from -1 to -15 listing number of +ve and -ve graphs
    with a weight in that band.  The top and bottom bands catch any overspills.


