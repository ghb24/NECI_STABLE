.. _input_non-block:

-----------------------
Non-block level options
-----------------------

The following options exist outside of any input block:

**TITLE**
   Takes the rest of the line as the title and prints it to output.  Useful for labelling the output.

**DEFAULTS** [ **DEFAULT** **FEB08** ]
   Default: **DEFAULT**.
   NECI has a default set of defaults (the **DEFAULT** set), which are sensible, safe defaults.  
   The **FEB08** set of defaults reflect furthr work, and change the defaults as follows:

       * Fock-Partition-Lowdiag is set in the integral block.
       * RhoEpsilon= :math:`10^{-8}` in the calc block.
       * MCPATHS is set to be on in the logging block. 

This can be specified anywhere in the input file outside of an input block.  All other options in the input file override the defaults.

**END**
   End of input file.  Not required, unless there is text after the input (e.g. comments or notes) which is not commented out or if the input file is given via STDIN.

