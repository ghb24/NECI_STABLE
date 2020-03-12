#!/usr/bin/env python

# this function histograms H_ij/p(i|j) to analyze the distribution of this ratio
# optimally this should have a maximum around 1 which means an optimal matrix element to 
# generation probability ratio. 
# input: 
# the input should be one or multiple files containing the matrix element and pgen information 
# and a histogram width if wished

# todo! 

from __future__ import (absolute_import, division, unicode_literals)
import sys
import glob
import re
import numpy
import matplotlib.pyplot as plt
from matplotlib import colors
from itertools import cycle

import six

if len(sys.argv) == 1: 
    quit("no data file provided! abort!")


# # if i want to analyze special excitation types:
# if len(sys.argv) > 2: 
#     spec_types = sys.argv[2:]
# 
#     spec_type = set([int(x) for x in spec_types])
# 
# else: 
#     # else use all types of excitations
#     spec_types = set(range(1,23))
# 
# 
# # make a list of markers for the different types of guga excitation if necessary
# # maybe also think about indicating similar excitation types with same marker or color.. TBD
mark = ['.',',','o','v','^','<','>','1','2','3','4','8','p','*','h','H','+','x','D','d']
# 
# # also create a color list
col = ['b','g','r','c','m','y','k']
# 
# # make a list of the specific names of excitations, to use them in the legend of the plot
typ_names = ['single', 'WR', 'WL','non-overlap','_L -> ^LL_ -> ^L', '_R -> ^RR_ -> ^R', '_L -> ^LR_ -> ^R', \
        '_R -> ^RL_ -> ^L', '_L -> _LL -> ^LL -> ^L', '_R -> _RR -> ^RR -> ^R', '_L -> _RL -> ^RL -> ^L', \
        '_R -> _LR -> ^LR -> ^R', '_L -> _RL -> ^LR -> ^R', '_R -> _LR -> ^RL -> ^L', '_L -> _LL -> ^LL^', \
        '_R -> _RR -> ^RR^', '_L -> _RL -> ^RL^', '_R -> _LR -> ^RL^', '_LL_ -> ^LL -> ^L', '_RR_ -> ^RR -> ^R', \
        '_RL_ -> ^LR -> ^R', '_RL_ -> ^RL -> ^L', '_RR_ -> ^RR^ / _LL_ -> ^LL^', '_RL_ -> ^RL^']

# i probably want to do that for all of the created pgen files or atleast have the option.. 

spec_types = set(range(1,23))
pgen = list()
matele = list()
typ = list()
# print sys.argv
for file_name in sys.argv[1:]:

    # file_name = sys.argv[1]
    tmp_file = open(file_name, 'r')

    # skip first line
    tmp_file.readline()

    csf = tmp_file.readline()

    # info about the processed CSF
    csf= re.findall(r'\d+', csf[0:csf.find(')')])

    # and convert that to integers
    csf = [int(x) for x in csf]

    # skip once more to get to data
    tmp_file.readline() 

    # now i have to determine the number of elements per row:
    # 2 ... no guga excitation type provided
    # 3 ... col 3 is guga information but no cummulative pgen info provided
    # 4 ... col 4 is guga information and cumulative pgen info is also provided 

    # do that by checking the first line
    line = tmp_file.readline()

    line = line.split()

    control = len(line)
# 
#     if control == 2:
# #         pgen = [float(line[0])]
# #         matele = [float(line[1])]
# 
#         for line in tmp_file:
#             line = line.split()
#             pgen.append(float(line[0]))
#             matele.append(float(line[1]))
# 
#         # convert them to numpy arrays
#         pgen = numpy.asarray(pgen)
#         matele = numpy.asarray(matele)
# 
#     elif control == 3: 
#         pgen = [float(line[0])]
#         matele = [float(line[1])]
#         typ = [int(line[2])]
# 
#         for line in tmp_file:
#             line = line.split()
#             pgen.append(float(line[0]))
#             matele.append(float(line[1]))
#             typ.append(int(line[2]))
# 
#     elif control == 4:
#     pgen = [float(line[0])]
#     matele = [float(line[1])]
#     typ = [int(line[3])]

    for line in tmp_file:
        line = line.split()
        pgen.append(float(line[0]))
        matele.append(float(line[1]))
        if (control >= 4): 
            typ.append(int(line[3]))

    tmp_file.close()

# no check the H_ij/pgen ratio and fill up bins 
# use numpy to do that 
# have to exclude all 0 from pgens.. 
matele = [matele[x] for x in range(len(matele)) if pgen[x] > 0.000001]

if control >= 4:
    typ = [typ[x] for x in range(len(typ)) if pgen[x] > 0.000001]

pgen = [x for x in pgen if x > 0.000001]

# pgen = [x for x in pgen if abs(x) > 0.0000001]
tmp_pgen = numpy.asarray(pgen)
tmp_mat = numpy.asarray(matele)

if control >= 4: 
    tmp_typ = numpy.asarray(typ)


hist, bin_edges = numpy.histogram(abs(tmp_mat)/tmp_pgen)
# print hist
# print bin_edges

plt.figure()
plt.hist(abs(tmp_mat)/tmp_pgen,50)
plt.draw()

# if i have excitation  typ information also make a histogram of the different type of excitations

if control == 4:

    # have to find all the different types of excitations
    types = set(typ)

    ind_list = range(0,len(pgen))
    # now plot the specific excitations with different markers and colors
    cnt = 0

#     min_rat = []
#     mean_rat = []

    # check if only certain excitations are of interest:
    for x in types & spec_types:
        plt.figure()
        tmp_pgen = [pgen[i] for i in ind_list if typ[i] == x and typ[i] in spec_types]
        tmp_mat = [abs(matele[i]) for i in ind_list if typ[i] == x and typ[i] in spec_types]

        # convert to numpy arrays:
        tmp_pgen = numpy.asarray(tmp_pgen)
        tmp_mat = numpy.asarray(tmp_mat)

        # also get the worst case ratio of pgen/hij to see, which type of excitations have the most detrimental 
        # effect on the time step 
#         min_rat.append(min(tmp_pgen/tmp_mat))

#         mean_rat.append(numpy.mean(tmp_pgen/tmp_mat))

        plt.hist(abs(tmp_mat)/tmp_pgen,50, color = col[cnt%7], label = str(x))
        plt.title(typ_names[x])
#         plt.scatter(abs(tmp_mat), tmp_pgen, c = col[cnt%7], marker = mark[cnt])
#         plt.grid()
#         plt.title("csf: " + str(csf))
        plt.legend()

        plt.draw()
        cnt = cnt + 1

plt.show()

