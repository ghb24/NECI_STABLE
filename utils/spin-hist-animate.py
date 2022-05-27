#!/usr/bin/env python

from __future__ import print_function
import sys, os
from subprocess import *

if __name__ == '__main__':

	if (len(sys.argv) < 4):
		sys.exit('Insufficient arguments given')

	# File name increment, and argument base
	fn_base = sys.argv[1]
	fn_inc = int(sys.argv[2])
	out_file = sys.argv[3]

	# Generate a sorted list of all the png files
	file_list = []
	iter_no = fn_inc
	fn = "%s%d.png" % (fn_base, iter_no)
	while os.path.isfile(fn):
		file_list.append(fn)
		iter_no += fn_inc
		fn = "%s%d.png" % (fn_base, iter_no)

	print("Collating %d png files" % len(file_list))

	Popen (["convert"] + file_list + [out_file]).communicate()
	






