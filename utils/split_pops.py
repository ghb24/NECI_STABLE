#!/usr/bin/python
'''Split and combine binary popsfiles for use with SPLIT-POPS

--> NEED the POPSFILEHEAD to be present for splitting the popsfiles up, as
    this is no longer such a 'dumb' process...

Usage:
    split_pops.py (combine [confirm]|split [num])

        combine - Combines the (consecutive) of the form POPSFILEBIN-[0-9]+
        confirm - Confirms overwriting an existing POPSFILEBIN
        split   - Splits a POPSFILEBIN into 'num' different ones for 'num' processors.
                  Will NOT overwrite existing POPSFILEBIN-0.'''

import struct
import sys
import re


def usage():
    '''Print the usage statement'''
    print __doc__




class fort_readwrite:

	def __init__ (self, f):

		# Pass in a file-like-object --> a bit more general than using
		# open directly here.
		self.f = f

		# Keep a rolling buffer of bits-and-pieces.
		self.buf = ""



	def read_to_buf (self):

		# Read the next chunk of the file into the buffer.
		tmp = self.f.read(4)
		if not tmp:
			return None
		(length,) = struct.unpack_from("@I", tmp)

		# Now read the actual data into the buffer
		self.buf += self.f.read(length)

		# Check that the 'end' block is correct.
		tmp = self.f.read(4)
		(length2,) = struct.unpack_from("@I", tmp)
		assert length == length2
		return length



	def write (self, buf):

		# We have to write the length in blocks at the start and end.
		length = len(buf)
		print "buflen", length
		tmp = struct.pack("@I", length)
		print 'tmplen', len(tmp)

		# Write the header, then the data, then the footer
		self.f.write(tmp)
		self.f.write(buf)
		self.f.write(tmp)



	def read(self, size):

		# Read until we have enough data in the buffer.
		while len(self.buf) < size:
			# If we have read all of the data, stop here.
			if not self.read_to_buf():
				break

		if len(self.buf) == 0:
			return None
		else:
			# Remove the returned data from the buffer.
			end = min(size, len(self.buf))
			ret = self.buf[0:end]
			self.buf = self.buf[end:]
			return ret


	def __enter__ (self):
		return self	



	def close(self):

		# Clear up the file that belongs to us.
		if self.f:
			self.f.close()
			self.f = None



	def __exit__ (self, type, value, traceback):

		print 'type', type
		print 'val', value
		print 'traceback', traceback
		print self.f
		self.close()


	def __del__ (self):

		self.close()




def ilut_to_ni (ilut, bits):
	'''Convet ilut to nI.'''

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



def signed_int_overflow (val, bits):
	'''To implement hashes, we need to emulate fortrans integer overflows'''


	# For hash calculation, we ensure use of 64-bit integers in NECI
	maxint = (1 << (bits - 1)) -1

	# We want to keep ourselves within the range of a _signed_ 64-bit integer
	if not -maxint-1 <= val <= maxint:
		val = (val + (maxint + 1)) % (2 * (maxint + 1)) - maxint - 1

	return val




def fmod (val, base):
	'''A modulus function that behaves like the fortran intrinsic'''

	tmp = val % base
	if val < 0 and tmp != 0:
		tmp = tmp - base

	return tmp




def determine_det_node (nI, random_hash, nnodes, bits):
	'''Determine which node determinant nI belongs to'''

	acc = 0
	for (i, orb) in enumerate(nI):
		tmp = (1099511628211 * acc) + (random_hash[orb - 1] * (i + 1))
		acc = signed_int_overflow(tmp, 64)

	return abs(fmod(acc, nnodes))





def process_header ():
	'''Process the header file. this will allow us to determine how many integers are used in each line...'''

	try:
		with open('POPSFILEHEAD', 'r') as f:

			nline = 0
			bits = 32
			random_hash = []
			for line in f:
			 
				txt = line.rstrip(', &\n').lstrip(', &\n')
				split = re.split(r'[, =]+', txt)
				nline += 1

				if nline == 1:

					# We are only really interested in popsfile version 4.
					if txt != "# POPSFILE VERSION 4":
						print "Invalid popsfile header"

				else:

					pos = 0
					while pos < len(split):

						if split[pos] == 'Pop64Bit':
							pos += 1
							bits = 64 if split[pos] == "T" else 32
							print "Using %d-bit POPSFILE" % bits

						elif split[pos] == "PopNIfTot":
							pos += 1
							niftot = int(split[pos])
							print "%d integers per determinant" % (niftot + 1)

						elif split[pos] == "PopNIfD":
							pos += 1
							nifd = int(split[pos])
							print "%d integers in spin-orb representation" % \
									            (nifd + 1)

						elif split[pos] == "PopRandomHash":
							pos += 1
							while pos < len(split):
								random_hash.append(int(split[pos]))
								pos += 1


						elif split[pos] == "END":
							print "End of header"
							break

						# Move on to the next entry
						pos += 1

	# If the popsfile header is not found, then there isn't much we can do.
	except IOError:
		print "POPSFILEHEAD not found"
		print ""
		usage ()
		sys.exit(-1)

	return (bits, niftot, nifd, random_hash)




def split_pops (nprocs):
	'''Split the popsfile into nprocs files of the format POPSFILEBIN-[0-9]+'''

	print "Splitting up POPSFILE bin into %d parts" % nprocs

	# Ensure that target POPSFILEBIN-* do not already exist
	for proc in range(nprocs):
		try:
			with open("POPSFILEBIN-%d" % proc, 'rb') as f:
				print "POPSFILEBIN-%d already exists" % proc
				print "Overwriting not supported."
				return
		except:
			pass

	# Extract the header information.
	(bits, niftot, nifd, random_hash) = process_header()

	# Open all of the output files
	outfiles = []
	for proc in range(nprocs):
		try:
			f = fort_readwrite(open("POPSFILEBIN-%d" % proc, 'wb'))
			outfiles.append(f)
		except:
			print "Error opening file: POPSFILEBIN-%d" % proc
			for f in outfiles:
				f.close()
			sys.exit(-1)

	# Read length (in bytes)
	readlen = (1 + niftot) * (bits / 8)

	# Loop through, reading 
	try:
		with fort_readwrite(open("POPSFILEBIN", 'rb')) as f:

			print "Opened"
			print "Read length: %d" % readlen

			unpack_str = "@" + ((nifd + 1) * "Q" if bits == 64 else "L") + "q"
			print "Unpack str |%s|" % unpack_str

			# Use a nifty sentinel value for iter, to read until there
			# is no more data!
			for cnt, data in enumerate(iter(lambda: f.read(readlen), None)):

				# Check that we have the right data length
				if len(data) != readlen:
					print "Mismatch in data lengths"
					sys.exit(-1)

				# Decode the determinant
				ilut_s = struct.unpack_from(unpack_str, data)
				ilut = ilut_s[0:nifd+1]
				sgn = ilut_s[-1]
				node = determine_det_node(ilut_to_ni(ilut, bits), random_hash, nprocs, bits)

				print ilut_to_ni(ilut, bits), node, sgn

				# Write the determinant back out to the relevant output file.
				outfiles[node].write(data)
				totwalkers = cnt + 1

			print "TOTWALKERS", totwalkers


			


	except IOError:
		print "Unable to open POPSFILEBIN"

	# Close all of the output files
	for f in outfiles:
		f.close()






def combine_pops (overwrite):
    '''Combine the popsfiles of the format POPSFILEBIN-[0-9]+'''

	# Don't automatically overwrite an existing POPSFILEBIN without
	# explicit confirmation.
    try:
        with open("POPSFILEBIN", 'rb') as fin:
            print "POPSFILEBIN already exists"
            if overwrite:
                print "OVERWRITING"
            else:
                print "Aborting"
                usage()
                return
    except:
        print "No existing POPSFILEBIN found. Creating it."


    # Open ourselves an output file
    with open("POPSFILEBIN", 'wb') as fout:

        # Loop through all of the present (and consecutive) POPSFILEBINs
        n = 0
        while True:

            fn = 'POPSFILEBIN-%d' % n
            try: 
                with open(fn, 'rb') as fin:
                    print "Opened %s" % fn

                    # Split the reads up into small chunks for memory usage
                    # control
                    while True:
                        data = fin.read(65536)
                        if not data:
                            break

                        # Write the data straight back to the combined FILE
                        fout.write(data)
                    
            except IOError:
                print "No more files"
                break

            n += 1




if __name__ == '__main__':

    
    # Determine what the input arguments are, and do the appropriate thing.
    if len(sys.argv) > 1 and sys.argv[1] == "combine":

        if len(sys.argv) > 2 and sys.argv[2] == "confirm":
            combine_pops(True)
        else:
            combine_pops(False)

    elif len(sys.argv) > 2 and sys.argv[1] == "split":
		split_pops (int(sys.argv[2]))
    else:
        usage()



