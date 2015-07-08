.. _output_mcpaths:

---------------------
MCPATHS and MCSUMMARY
---------------------

The MCPATHS logging file (or MCSUMMARY if a **METHODS** block is used) has the following layout::

    <Header Line>
    <Vertex Sum Section for Det 1>
    <Vertex Sum Section for Det 2>
    ...
    <MC Summary information>

The Header Line is::

    \"Calculating  XXX W_Is...\"

where XXX is the number of determints calculated, and the number of Vertex Sum sections.

Each Vertex Sum Section consists of::

        (Determinant Calculated)
        <method level 1 line>
        <method level 2 line>
        ...

The Method level line is::

    <level>  <weight>  <cumlweight>  <timing> <GraphsSummed>  [<PartGraphs>] <w E~> [<MP2 contrib> [<MP3 contrib> [..]]]]

where:
    level
        The vertex level as specified by the METHODS section.
    weight
        The contribution of this level to s_i (see RHOPII or RHOPIIex file section).  This can be further analysed in CLASSPATHS{,2}.
    cumlweight
        The sum of weights up to and including this level.
    timing
        (double) The number of seconds calculating this level took.
    GraphsSummed
        The total number of graphs summed together at this level.
    PartGraphs
        Recursively, how many times FMCPR?  is called.  It is called once as each node is added to a graph.
    :math:`w \tilde{E}`
        The contribution of this level to :math:`w \tilde{E}`.  This can be further analysed in CLASSPATHS{,2}
    MPn contrib
        The contribution of graphs at this level to MPn theory.

There are is one method level line for each method level specified in the **METHODS** section, plus one  for the 1-vertex graph.

The MC Summary information is split into 4 parts:

    TotalStats::

        GRAPHS(V)...WGHT-(V)

    GenStats::

        GEN-> *

    AccStats::

        ACC-> *

    Sequences::

        Sequences, Seq Len

TotalStats:
    The output is split into columns depending on the levels sampled in the Monte Carlo:

    DataType    Total    1-vertex    2-vertex    3-vertex    ...

    The DataTypes are:
        GRAPHS(V)
		    The number of graphs sampled with this number of vertices.
        TREES(V)
		    The number of trees sampled with this number of vertices.  Trees contain no cycles.
        NON-TR+(V)
		    The number of non-trees sampled with this number of vertices whose weight is positive.
        NON-TR-(V)
		    The number of non-trees sampled with this number of vertices whose weight is negative.
        WGHTT(V)
		    The total weight of trees with this number of vertices.
        WGHT+(V)
		    The total weight of positive non-trees with this number of vertices.
        WGHT-(V)
		    The total weight of negative non-trees with this number of vertices.

GenStats:
    Statistics on the number of graph-graph transitions generated.

    The columns correspond to the graph FROM which the transition was generated:
        DataType    1-vertex    2-vertex    3-vertex    ...
    The rows correspond to the graph TO which the transition was generated:
        GEN-> 1
        GEN-> 2
        ...

AccStats:
    Format as GenStats, but has the number of transitions accepted.

Sequences:
    Records the number of sequences of consecutive graphs accepted with the same weight.

