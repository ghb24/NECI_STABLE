#!/usr/bin/env python

from __future__ import print_function
import sys

#
# N.B. ALL of the following assume unsigned integers
#
def ilut_to_ni (ilut, bits):
	'''Convert ilut to nI. Assume ilut is length 1 for now.'''

	nI = []
	offset = 0
	for tmp in ilut:
		i = 0
		while (1 << i) < tmp:
			if not (tmp & (1 << i)) == 0:
				nI.append(i + 1 + offset)
			i += 1
		offset += bits

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

def count_bits (ilut):
	'''Naive bit counting algorithm'''

	n = 0
	for tmp in ilut:
		i = 0
		while (1 << i) < tmp:
			if (tmp & (1 << i)) != 0:
				n += 1
			i += 1

	return n



class sorted_list:
	def __init__ (self, nmax):
		self.keys = []
		self.values = []
		self.nmax = nmax

	def insert (self, key, value):
		'''Insert a value into the list ONLY if big enough'''
		if len(self.values) < self.nmax or value > self.values[0]:

			# If we need to shrink the list appropriately
			if len(self.values) >= self.nmax:
				del self.values[0]
				del self.keys[0]

			for i in range(len(self.values)):
				if self.values[i] > value:
					self.values.insert(i, value)
					self.keys.insert(i, key)
					return

			self.values.append(value)
			self.keys.append(key)

	def sorted_list (self):
		ind = 0
		for key in self.keys:
			yield (key, self.values[ind])
			ind += 1




# What popsfile should we use?
if __name__ == "__main__":
	if len(sys.argv) > 1:
		pops = sys.argv[1]
	else:
		pops = 'POPSFILE'
	print('Using POPS file: %s' % pops)

	# Parameters
	nfind = 20
	nopen_req = 0
	update_iters = 10000

	# Dictionary of largest values
	largest = sorted_list(nfind)
	nlargest = 0

	nw_found = 0
	with open(pops, 'r') as f:

		# Read all the lines in the file
		nline = 0
		for line in f:

			txt = line.rstrip()
			nline += 1

			# Check version/header information
			if nline == 1:
				if txt == "# POPSFILE VERSION 3":
					ver = 3
				elif txt == "# POPSFILE VERSION 2":
					ver = 2
				else:
					sys.exit('Invalid popsfile version')
			elif ver == 3:
				if nline == 2:
					bits = 64 if bool(txt.split()[1]) else 32
					nel = int(txt.split()[9])
				elif nline == 3:
					nw = int(txt.split()[0])
				elif nline == 8:
					nifd = int(txt.split()[0])
				elif nline == 9:
					nify = int(txt.split()[0])
				elif nline == 10:
					nifsgn = int(txt.split()[0])
				elif nline == 11:
					nifflg = int(txt.split()[0])
				elif nline == 12:
					niftot = int(txt.split()[0])
					break

			elif ver == 2:
				niftot = -1
				nify = -1
				nifflg = -1
				nifsgn = -1
				nifd = -1
				nel = -1
				if nline == 2:
					bits = 64 if bool(txt.split()[1]) else 32
				elif nline == 3:
					nw = int(float(txt.split()[0]))
				elif nline == 7:
					break

		# Output the gathered information
		# (These should all be defined by now...)
		print('nel: %d' % nel)
		print('bits', bits)
		print('nwalkers: %d' % nw)
		print('NIfD: %d' % nifd)
		print('NIfY: %d' % nify)
		print('NIfSgn: %d' % nifsgn)
		print('NIfFlag: %d' % nifflg)
		print('NIfTot: %d' % niftot)
		print('----------------')

		# unsigned value of -1:
		minus1 = 0xFFFFFFFF if bits == 32 else 0xFFFFFFFFFFFFFFFF

		for line in f:

			nline += 1

			if nw_found % update_iters == 0:
				print('%d/%d: %d' % (nw_found, nw, nlargest))

			# Is there an incomplete line at the end of the POPSFILE
			try:
				tsplit = map(int, line.rstrip().split())

				# Make a stab at this information if we need to
				if ver == 2 and niftot == -1:
					nifsgn = 1
					nifflg = 0
					nify = 0
					niftot = len(tsplit) - 1
					nifd = niftot - nifsgn
					wtmp = map(lambda x: x if x >= 0 else minus1+x+1, 
					           tsplit[0:nifd+1])
					nel = count_bits(wtmp)
					print('UPDATED')
					print('nel: %d' % nel)
					print('bits', bits)
					print('nwalkers: %d' % nw)
					print('NIfD: %d' % nifd)
					print('NIfY: %d' % nify)
					print('NIfSgn: %d' % nifsgn)
					print('NIfFlag: %d' % nifflg)
					print('NIfTot: %d' % niftot)
					print('----------------')

			except:
				print('Invalid walker found in popsfile, line %d' % nline)
				sys.exit('bad')
				continue

			if len(tsplit) < niftot + 1:
				print('Incomplete walker found in popsfile, line %d' % nline)
				continue

			# Ensure we use unsigned integers for these.
			w = map(lambda x: x if x >= 0 else minus1+x+1, tsplit[0:nifd+1])
			sgn = tsplit[nifd+1]
			nlargest = max(nlargest, abs(sgn))
			nop = nopen(w)

			# Keep track of number of walkers found
			nw_found += 1

			# Should we store this one?
			if nop >= nopen_req:
				largest.insert(w, abs(sgn))

		print('----------------')
		print('nlargest', nlargest)
		print('Largest items')
		for pair in largest.sorted_list():
			print(pair[0], end=' ')
			print('%d - ' % pair[1], end=' ')
			print(ilut_to_ni(pair[0], bits))

		print('Cleaning up')
		f.close()

