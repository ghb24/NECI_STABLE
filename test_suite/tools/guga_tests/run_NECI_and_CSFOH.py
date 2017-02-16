#!/usr/bin/env python

# this script should be able to run neci followed by a matrix element comparison with Block! 
# this will be used in the testsuite to provide a fully automatic testsuite for my GUGA branch, which checks 
# that the pgens and the matrix elements are correct!

import subprocess
import os
from numpy.testing import assert_allclose
import numpy
# from comparing_mod import *
import glob
import re
import sys
import socket

# first of all clean the directories: especially from any created contribs_guga.* and pgen_* files
map(os.remove, glob.glob('contribs_guga.*'))
map(os.remove, glob.glob('pgen_vs_*'))

# try to run a neci calc and fetch the output:
# i have to ensure that the neci process finishes! so set a nmcyc! and i have to deal with the case that the 
# testsuite does a stop_all()! test that! 
# make this script now more general to run on any machine i am using
machine = socket.gethostname()
if machine == 'baloo-X1-Carbon':
    neci_exe = '/home/dobrautz/bin/neci.x'
    compare_exe = '/home/dobrautz/cloud/doktorat/neci/neci/test_suite/tools/guga_tests/compare_all_excitations.py'

elif machine == 'pcal008': 
    neci_exe = '/home/dobrautz/bin/neci.x'
    compare_exe = '/home/dobrautz/cloud/doktorat/neci/neci/test_suite/tools/guga_tests/compare_all_excitations.py'

elif machine == 'allogin1':
    neci_exe = '/home/dobrautz/bin/neci.cluster.x'
    compare_exe = '/home/dobrautz/cloud/doktorat/neci/neci/test_suite/tools/guga_tests/compare_all_excitations.py'

elif machine == 'altest':
    neci_exe = '/home/dobrautz/bin/neci.x'
    compare_exe = '/home/dobrautz/cloud/doktorat/neci/neci/test_suite/tools/guga_tests/compare_all_excitations.py'

else:
    print 'incorrect hostname! check setup!'
    raise SystemExit()

# check if exe exists
if not os.path.isfile(neci_exe): 
    print 'neci exe not found!'
    raise SystemExit()

if not os.path.isfile(compare_exe): 
    print 'neci exe not found!'
    raise SystemExit()
process_neci = subprocess.Popen([neci_exe,'input.inp'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

output_neci, error_neci = process_neci.communicate()

print output_neci
print error_neci

# yes that works! now i have to check how to get the comparison running depending if the neci calc. was 
# successful!
# do some testing if neci did not throw an error!
# ok neci does not quite properly work with the errorcodes.. so i could try to check if the error_message is empty:
if error_neci == '': 
    # if neci worked do the dmrg_comparison
    print "******************************************"
    print " finished NECI calculation!"
    print "NECI Run successful:", 1
    print " Now entering CSFOH matrix element comparison"
    print "******************************************"

    # and then try to call my compare_all python script.. 

    process_dmrg = subprocess.Popen([compare_exe], stdout=subprocess.PIPE, \
            stderr=subprocess.PIPE)

    output_dmrg, error_dmrg = process_dmrg.communicate()

    # i have to change the output off my compare_all routine to give me something meaningful to compare with 
    # actually i can stop the comparison if the first wrong matrix element is found and exit with a negative 
    # exit code or smth.
    print output_dmrg
    print error_dmrg

    print "still running?"

    # i also have to check if CSFOH did not crash! due to bad FCIDUMPs or some stuff! 
    if error_dmrg == '':
        if process_dmrg.returncode == 1:
            # output that the tests failed, i have to indicate that neci worked, but dmrg did not! 
            print "*******************************************"
            print " CSFOH ran, but matrix elements are incorrect!"
            print "CSFOH Run succesful:", 0
            print "*******************************************"

        else:
            # the neci and the dmrg comparison worked! so combine the output that it can be sure that neci and dmrg 
            # worked! 
            print "*******************************************"
            print " CSFOH ran and matrix elements are correct!"
            print "CSFOH Run succesful:", 1
            print "*******************************************"

    else:
        # CSFOH crashed!

        print "*******************************************"
        print " CSFOH crashed with error message:"
        print error_dmrg
        print "CSFOH Run succesful:", 0
        print "*******************************************"

else: 
    # the neci calc failed already: indicate that in the output, so that testcode.py picks it up and gives a failed 
    # test result!
    print "*******************************************"
    print " NECI calculation aborted with error: "
    print error_neci
    print "NECI Run successful:", 0
    print " Cannot enter CSFOH matrix element comparison!"
    print "CSFOH Run succesful:", 0
    print "*******************************************"

# some cleaning up
map(os.remove, glob.glob('contribs_guga.*'))
map(os.remove, glob.glob('pgen_vs_*'))
