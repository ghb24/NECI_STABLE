#!/bin/sh

# little script to run a minimal test_suite in bitbucket pipelines to save some build time
# the default installation path of testcode in my docker image of the neci setup: dobrautz/neci_base 
# is /testcode/bin/testcode.py 

/testcode/bin/testcode.py -c dneci/double_occ/hub_2x2 -c kmneci/Rn_lanczos_dci_init -c kneci/C_221_int -c mneci/cfqmc/HeHe_5_states -c neci/parallel/C_Solid
