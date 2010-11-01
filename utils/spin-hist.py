#!/usr/bin/python

'''Plot a histogram from a spin-hist-* file. This currently assumes that the system is only using real walkers.'''

import sys
import pylab
import os, re
from subprocess import *
from math import *

if __name__ == '__main__':

	if (len(sys.argv) < 3): #4):
		sys.exit('Insufficient arguments given')

	# File name increment, and argument base.
	fn_base = sys.argv[1]
	fn_inc = int(sys.argv[2])
	#stats_file = sys.argv[3]

	# What stats files are there?
	test = re.compile("%s[0-9]*$" % fn_base)
	command = filter(test.match, os.listdir(os.getcwd()))
	command = ["awk", "($2 > 0 ? $2 : 0-$2) > max && NR /= 1 {max = ($2 > 0 ? $2 : 0-$2)} END{print max}"] + command

	# Get the maximum count that we need
	max_val = Popen(command, stdout=PIPE).communicate()[0].strip()
	max_val = int(ceil(float(max_val)))
	print 'Absolute max value: %d' % max_val

#	# Process the stats file to get totwalkers
#	tot_walkers = [0]
#	with open(stats_file, 'r') as f:
#
#		# Get iteration and walker number counts.
#		it, wlk = pylab.loadtxt (f, usecols=(0, 4), unpack=True)
#
#		itr = fn_inc
#		for i in range(len(it)):
#
#			if it[i] >= itr:
#				tot_walkers.append(wlk[i])
#				max_wlk = wlk[i]
#				itr += fn_inc


	mult = 1
	cum_max = 1.0
	while True:
		fn = "%s%d" % (fn_base, mult * fn_inc)

		# Leave the loop when we are done
		if not os.path.isfile(fn):
			break

		print 'Iter: %d (max: %f)' % (mult * fn_inc, cum_max)
		with open(fn, 'r') as f:
			sgn = []
			firstline = True
			for ln in f.readlines():

				if firstline:
					firstline = False
					continue

				ltmp = ln.split()
				if ltmp[0][0] == '#':
					continue

				#sgn.append(max_wlk * float(ltmp[len(ltmp)-1]) / tot_walkers[mult])
				sgn.append(float(ltmp[len(ltmp)-1]))
			
			# Maintain cumulative maximum
			cum_max = max(cum_max, max(map(abs, sgn)))

			fig = pylab.figure()
			ax1 = fig.add_subplot(212)
			ax2 = fig.add_subplot(211)

			ax1.bar(range(len(sgn)), sgn, width=1.0)
			ax1.set_xlim(0, len(sgn))
			ax2.bar(range(len(sgn)), sgn, width=1.0)
			ax2.set_xlim(0, len(sgn))
			ax1.set_ylim(-max_val, max_val)
			ax2.set_ylim(-cum_max, cum_max)
			pylab.savefig('%s%s' % (fn, '.png'))
#			pylab.show()
#			pylab.clf()

		mult += 1





