#!/usr/bin/python

import sys

def ilut_to_ni (ilut):
	'''Convert ilut to nI. Assume ilut is length 1 for now.'''

	i = 0
	nI = []
	while (1 << i) < ilut:
		if not (ilut & (1 << i)) == 0:
			nI.append(i + 1)
		i += 1
	return nI

def nopen (ilut):
	'''Determine the number of unpaired electrons'''

	n = 0
	for tmp in ilut:
		i = 0
		while (1 << i) < tmp:
			if (((tmp & (1 << i)) == 0) is not ((tmp & (1 << i+1)) == 0)):
				n += 1
			i += 2
	return n


# What popsfile should we use?
if __name__ == "__main__":
	if len(sys.argv) > 1:
		pops = sys.argv[1]
	else:
		pops = 'POPSFILE'
	print 'Using POPS file: %s' % pops

	nfind = 20
	nopen_req = 4
	with open(pops, 'r') as f:

		# Check version information
		for line in range(12):
			ltmp = f.readline().rstrip()

			if line == 0:
				assert ltmp == "# POPSFILE VERSION 3"
			elif line == 1:
				bits = 64 if bool(ltmp.split()[1]) else 32
				nel = int(ltmp.split()[9])
			elif line == 2:
				nw = int(ltmp.split()[0])
			elif line == 7:
				nifd = int(ltmp.split()[0])
			elif line == 8:
				nify = int(ltmp.split()[0])
			elif line == 9:
				nifsgn = int(ltmp.split()[0])
			elif line == 10:
				nifflg = int(ltmp.split()[0])
			elif line == 11:
				niftot = int(ltmp.split()[0])

		# Output the gathered information
		print 'nel: %d' % nel
		print 'bits', bits
		print 'nwalkers: %d' % nw
		print 'NIfD: %d' % nifd
		print 'NIfY: %d' % nify
		print 'NIfSgn: %d' % nifsgn
		print 'NIfFlag: %d' % nifflg
		print 'NIfTot: %d' % niftot


		largest = {}
		top = 0
		for i in range(nw):
			
			if i % 100 == 0:
				print '%d: %d' % (i, len(largest))

			tmp = f.readline().rstrip().split()

			w = map(int, tmp[0:nifd+1])
			sgn = int(tmp[nifsgn])
			top = max(top, abs(sgn))
			nop = nopen(w)

			if nop >= nopen_req and (len(largest) == 0 or abs(sgn) > min(largest.values())):

				largest[w[0]] = abs(sgn)
				if len(largest) > 10:
					key = None
					valmin = None
					for k,v in largest.items():
						if not valmin or v < valmin:
							valmin = v
							key = k
					del largest[key]


		print 'Largest items'
		sortkeys = largest.keys()
		sortkeys.sort()
		for k in sortkeys:
			print '(',
			print  k,
			print ') %d - ' % largest[k],
			print ilut_to_ni(k)

		print 'top', top

		print 'Cleaning up'
		f.close()

