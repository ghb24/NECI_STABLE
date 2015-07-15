.. _output_classpaths:

----------
CLASSPATHS
----------

CLASSPATHS is calculated when a vertex sum is performed and lists properties of the graphs by their class::

    Class nGs TotWeight   TotwEt TotWeightPos   TotWeightNeg

where:

    Class is a binary string indicating the connectivity of the graph.

    The connectivity matrix has 1 if there's a line in the graph
       e.g. a 3-vertex graph.  The bits are labelled from the bottom right starting from 0::
 
          (. 2 1 )  (bits)                   (. 1 1 )
          (  . 0 )  -> 210     connections:  (  . 0 ) -> 110 (binary) -> 6 (decimal)
          (    . )                           (    . ) 

    nGs is the number of graphs that fell in this class.
    TotWeight is the total weight, :math:`w_i`, of those graphs.
    TotEt is the total value of :math:`w_i \tilde{E}_i` for those graphs.
    TotWeightPos is the total weight of graphs(?) with positive weights.
    TotWeightNeg is the total weight of graphs(?) with negative weights.


