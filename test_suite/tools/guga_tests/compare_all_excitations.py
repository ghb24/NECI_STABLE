#!/usr/bin/env python

# script to compare matrix elements created by a NECI-GUGA run and Sandeeps Block code

import subprocess
import os
from numpy.testing import assert_allclose
import numpy
from comparing_mod import *
import glob
import re
import sys
import socket

# make this file more flexible for different systems
machine = socket.gethostname()
if machine == 'baloo-X1-Carbon':
    csfoh_exe = '/home/dobrautz/bin/CSFOH'

elif machine.startswith('pcal'):
    csfoh_exe = '/home/dobrautz/bin/CSFOH'

elif machine == 'allogin1':
    csfoh_exe = '/home/dobrautz/bin/CSFOH.cluster'

elif machine == 'altest':
    csfoh_exe = '/home/dobrautz/bin/CSFOH'

else:
    print 'incorrect hostname! check setup!'
    raise SystemExit()

# check if exe exists
if not os.path.isfile(csfoh_exe): 
    print 'CSFOH exe not found!'
    raise SystemExit()


# have to read in a output file of the guga_testsuite -> use the name of the file as an input 
# or loop over all possible files in the current folder and do the test on all of them

# from those files read the CSFs, store the matrix elements and produce input files for 
# CSFOH to run on

# and then run CSFOH and compare the matrix elements 


# if no input is provided check the current directory for all files named contribs_guga.*
if len(sys.argv) == 1 or sys.argv[1] == '':
    fileList = glob.glob("contribs_guga.*")
else:
    fileList = sys.argv[1:]

if len(fileList) == 0:
    quit("no contribs_guga.* files found! Aborting!")

print "*******************************"
print " Comparing matrix elements for output files:"
print fileList

# then loop over all the found or provided files
for filename in fileList: 

    print 
    print " Processing: ", filename
    file = open(filename, 'r')

    # skip the firt line
    file.readline()

    csf = file.readline()

    # this gives a list of the original stepvectors 
    csf = re.findall(r'\d+', csf[0:csf.find(')')])

    # and convert to ints
    csf = [int(x) for x in csf]

    # some output
    print " Comparing GUGA and DMRG matrix elements for CSF: "
    print csf

    # skip one more line to get to the excitations
    file.readline()

    excitations = []
    mat_eles = []
    # loop over the excitations and extract csfs and matrix elements 
    for line in file: 
        # make a list of the excitations
        strings = re.findall(r'\d+', line[0:line.find(')')])
        excitations.append([int(x) for x in strings])

        # and of the corresponding matrix elements
        mat_eles.append(float(line.split()[-1]))

    # close the file
    file.close()

    print " GUGA Excitations: "
    for i, item in enumerate(excitations):
        print item, mat_eles[i]


    print " Writing DMRG config and determinants file"
    # and then write the DMRG config and determinants file to run CSFOH
    writeDMRGconfig(csf)


    ref_mat_eles = []
    cnt = 0
    # up until here it is fine i guess, now i just need to compare the first CSF with all the other ones one at a 
    # time! loop over them! so i want to do one DMRG calc for every excitation and the reference alone! 
    for excit in excitations:
        excit = [excit]
        writeDeterminants(csf, excit)

        # then run CSFHOF
        process = subprocess.Popen([csfoh_exe,'dmrg.conf'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

#         print " Running DMRG calculation ..."
        output, error = process.communicate()

#         print " Finished DMRG calculations!"
        # find correpoding weights
        start = output.find('>\n')
        end = output.find('\n',start + 3)

        ref_mat_ele = output[start + 3: end]

#         print output
        # convert to floats
#         print ref_mat_ele
        ref_mat_ele = [float(x) for x in ref_mat_ele.split()]
        # the second entry is the overlap matrix element! 
        ref_mat_eles.append(ref_mat_ele[1])
        cnt += 1
#         print "iteration: ", cnt
        

    print 
    print " Comparing matrix elements:"
    print " DMRG | GUGA "
    # compare the matrix elements 
    cnt = 0
    for i, item in enumerate(ref_mat_eles):
        print item, mat_eles[i]
        if abs((abs(item) - abs(mat_eles[i]))/abs(mat_eles[i])) > 0.001:
            print " Incorrect matrix element for excitation:"
            print excitations[i]
            print " DMRG: ", item, " GUGA: ", mat_eles[i]
            cnt += 1

    if cnt == 0: 
        print "======================================"
        print " All matrix elements correct for CSF: "
        print csf
        print "======================================"

    else:
        print "======================================"
        print " NOT all matrix elements correct!"
        print "======================================"


        # cleaning up
        print "cleaning up"
        os.remove('determinants')
        os.remove('dmrg.conf')

        # and exit with
        sys.exit(1)

# cleaning up
print "cleaning up"
os.remove('determinants')
os.remove('dmrg.conf')
sys.exit(0)


