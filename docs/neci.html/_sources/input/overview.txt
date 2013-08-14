.. _input_overview:

--------
Overview
--------

The NECI input file is keyword driven and requires a minimal amount of information.  

The NECI input file is divided into various sections, or input blocks: system, precalc, calc, integral and logging.  Of these, only the system and calc blocks are compulsory: all others are optional.  Inside each input block, it is possible to set a variety of options.  There are also three types of keywords that exist outside of an input block.

The order of the input blocks is not important (but certain orders are more logical than others), and nor is the order within a block, unless an option is only valid when a logical statement is true, in which case the relevant keyword for the logical statement must precede its related keywords.

General points to note:

* The input file is not case sensitive.  In the input documentation, the keywords are given in capitals and **emphasised** for clarity and options or data required are in square brackets.
* Parameters which follow a keyword ought to be on the same line as the keyword,but this isn't a strict requirement.
* A new line is required for each keyword, unless the keyword is an option of another keyword, in which case it ought to be on the same line.
* Blank lines are ignored.
* Comments are enclosed in parentheses.
* Data items are terminated by space or comma.
* Only the variables relevant to the desired run are required.
* Unknown keywords return an error message and stop the run.
* Sensible defaults are set, reducing the amount of information required from the input file.  There exist different sets of default options, allowing a large set of variables to be set with one command.

The overall structure, with a reasonably logical layout, is:

**TITLE**

**DEFAULTS**

**SYSTEM** [system type]

[System options]

**ENDSYS**

**PRECALC**

[PreCalc options]

**ENDPRECALC**

**CALC**

[Calc options]

**ENDCALC**

**INTEGRAL**

[Integral options]

**ENDINT**

**LOGGING**

[Integral options]

**ENDLOG**

**END**

.. warning::
  This is a work in progress.  Many places (especially, but not
  exclusively, where noted) need to be expanded and/or improved.

  In addition, the following keywords are valid options, but are
  *not* documented:

     *  CALCREALPROD
     *  CALCRHOPROD
     *  DELTAH
     *  DERIV
     *  DETPOPS
     *  DIAGSHIFT
     *  EQUILSTEPS
     *  EXCHANGE-ATTENUATE
     *  HAPP
     *  LINROOTCHANGE
     *  MAXVERTICES
     *  MODMPTHEORY
     *  RESUMFCIMC
     *  RHOAPP
     *  RHOELEMS
     *  SAVEPREVARLOGGING
     *  SHIFTDAMP
     *  STARPROD
     *  STEPSSHIFT
     *  SUMPRODII 

   In contrast, the following options are documented, but are *not* valid
   input options:

      *  BANDGAP
      *  EXCHANGE-DAMPING
      *  STOCHASTICTIME
      *  MPMODTHEORY
      *  SAVEPRECALCLOGGING 
