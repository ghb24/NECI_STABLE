#!/usr/bin/python

'''Plot a histogram from a spin-hist-* file. This currently assumes that the system is only using real walkers.

Usage: spin-hist.py [args]

Options:
	--help, -h          Display this message
	--hist-base, -b     File name base used. default: spin-hist-
	--iter-inc, -i      Iteration increment. Use iterations inc, 2*inc...
	                    default: 100
	--plot-csfs, -c     Plot csfs instead of absolute amplitude as lower
	                    plot. defalt: False
	--csf-base, -C      Base filename for csf plots. default: csf-hist-
'''

import pylab
import getopt
import os, re, sys

from subprocess import *
from math import *


def usage ():
	'''Print the usage statement'''
	print __doc__


if __name__ == '__main__':

	# Set defaults
	fn_base = 'spin-hist-'
	fn_inc = 100
	csf_base = 'csf-hist-'
	plot_csf = False

	# Process input arguments
	try:
		opts, unused = getopt.getopt(sys.argv[1:], "b:i:cC:", ["hist-base=", "iter-inc=", "plot-csfs", "csf-base="])
	except getopt.GetopetError, err:
		print str(err)
		usage()
		sys.exit(2)

	for o, a in opts:
		if o in ("-h", "--help"):
			usage()
			sys.exit()
		elif o in ("-b", "--hist-base"):
			fn_base = a
		elif o in ("-i", "--iter-inc"):
			fn_inc = int(a)
		elif o in ("-c", "--plot-csfs"):
			plot_csf = True
		elif o in ("-C", "--csf-base"):
			csf_base = a
		else:
			print 'Unhandled option: %s' % o
			usage()
			sys.exit(2)


	# We only need to find the absolute maximum value if we are plotting
	# an overall plot. We don't do this if we have csfs.
	if not plot_csf:

		# What histogram files are there?
		test = re.compile("%s[0-9]*$" % fn_base)
		command = filter(test.match, os.listdir(os.getcwd()))
		command = ["awk", "($2 > 0 ? $2 : 0-$2) > max && NR /= 1 {max = ($2 > 0 ? $2 : 0-$2)} END{print max}"] + command

		# Get the maximum count that we need
		max_val = Popen(command, stdout=PIPE).communicate()[0].strip()
		max_val = int(ceil(float(max_val)))
		print 'Absolute max value: %d' % max_val


	mult = 1
	cum_max = 1.0
	cum_max_csf = 1.0
	while True:
		fn = "%s%d" % (fn_base, mult * fn_inc)

		# Leave the loop when we are done
		if not os.path.isfile(fn):
			break

		print 'Iter: %d (max: %f)' % (mult * fn_inc, cum_max)
		with open(fn, 'r') as f:

			# Extract histogram data
			sgn = []
			firstline = True
			for ln in f.readlines():

				if firstline:
					firstline = False
					continue

				ltmp = ln.split()
				if ltmp[0][0] == '#':
					continue

				sgn.append(float(ltmp[len(ltmp)-1]))

			if plot_csf:
				fn_csf = "%s%d" % (csf_base, mult * fn_inc)
				with open(fn_csf, 'r') as fcsf:
					csf_sgn = pylab.loadtxt(fcsf, usecols=[0], unpack=True)
					cum_max_csf = max(cum_max_csf, max(map(abs, csf_sgn)))
			
			# Maintain cumulative maximum
			cum_max = max(cum_max, max(map(abs, sgn)))

			fig = pylab.figure()
			ax2 = fig.add_subplot(212)
			ax1 = fig.add_subplot(211)

			ax1.bar(range(len(sgn)), sgn, width=1.0)
			ax1.set_xlim(0, len(sgn))
			ax1.set_ylim(-cum_max, cum_max)

			if plot_csf:
				ax2.bar(range(len(csf_sgn)), csf_sgn, color='r', width=1.0)
				ax2.set_xlim(0, len(csf_sgn))
				ax2.set_ylim(-cum_max_csf, cum_max_csf)
			else:
				ax2.bar(range(len(sgn)), sgn, width=1.0)
				ax2.set_xlim(0, len(sgn))
				ax2.set_ylim(-max_val, max_val)

			pylab.savefig('%s%s' % (fn, '.png'))
#			pylab.show()
#			pylab.clf()

		mult += 1





