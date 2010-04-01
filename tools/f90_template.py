#!/usr/bin/python
'''Produce a .f90 file from the given .F90.template file to simulate the
effects of using templates/generic functions in a language such as c++

Usage:
	f90_template.py [infile] [outfile]

TODO: fill in format information
'''
import os
import sys
import ConfigParser
import re

class file_like_list(list):
	'''Implement a list with  the readline function, so that it can act as
	though it were a file object.'''

	def __init__(self):
		self.pos = 0

	def readline(self):
		self.pos += 1
		if self.pos > self.__len__():
			return None
		else:
			return self[self.pos-1]


def usage():
	'''Prints the usage description for the program'''
	print __doc__

def read_config(fin):
	'''Read in the config from the start of the template file. Name/value
	pairs are generally preserved between sections unless overridden.'''

	config = {}
	config_lines = file_like_list()
	sections = list()

	# Read in the config data to a file like object, so that ConfigParser can
	# work on less than a whole file.
	re_section = re.compile('^\[(.*)\]$')
	while True:
		line = fin.readline()
		if line[0] == '=':
			break
		config_lines.append(line)

		# Python 2.6 does not have an ordered dictionary. Therefore to
		# preserve section order, we have to read it in.
		m = re_section.match(line)
		if m:
			sections.append(m.group(1))


	# Parse the first section of the config file for templated settings to
	# insert into the templated functions.
	parser = ConfigParser.RawConfigParser()

	parser.readfp(config_lines)
	
	print 'Producing configurations: ',
	sprev = None
	for s in sections:
		if sprev is not None:
			config[s] = dict()
			config[s].update(config[sprev])
			config[s].update(parser.items(s))
		else:
			config[s] = dict(parser.items(s))
		config[s].update({'name' : s})

		sys.stdout.write("%s%s" % (", " if sprev is not None else "", s))
		sprev = s
	print

	return config

def super_module(template, config):
	'''Generate a super module, which simply includes all of the variations
	on the template module. This routine also changes the template so that the
	specified module name will be varied in the output (so that it can be
	properly included in the super-module).'''

	re_mod = re.compile ('(\n\s*module[\s^\n]*([^\s]*))\n')

	m = re_mod.search(template)
	print 'Generating super module: %s' % m.group(2)

	# Change module in template to be a sub-module with associated name
	substr = "%s_%%(name)s\n" % m.group(1)
	template = re_mod.sub (substr, template)

	# Construct super module
	super_mod = m.group(1) + "\n"
	for s in config:
		super_mod += "    use " + m.group(2) + "_" + s + "\n"
	super_mod += "end module\n"

	return (template,super_mod)

def adj_arrays (template, config):
	'''Find all of the types we are specifying, and adjust references to their
	arrays in the code so we can vary the dimensionality of a template fn.
	Only tested for dimensions 1,2 at the moment, so most likely breaks at
	higher dimensionality.'''

	# All possible types we are looking for. Only interested if the
	# dimensionality > 1 (i.e. we have to adjust the dimensions)
	types = dict()
	re_dim = re.compile ('([\s,]*)dimension(\(.*\))')
	for i in config:
		if i[0:4] == "type":
			if config[i][0] != "!":
				m = re_dim.search(config[i])
				if m:
					types[i] = m.group(2).count(":")
					config[i] = re_dim.sub("", config[i])
				else:
					types[i] = 0


	# What variables occur with those types? Store the dimensionality.
	# I am not happy with how general this is - I am sure there are things
	# which will break it.
	#     --> TODO: reconsider
	#     --> e.g. If we declare multiple variables to be reshaped on the
	#              same line, then it wont work.
	#     --> TODO: restrict the scope of each search and replace to one fn.
	
	# Vars contains tuples of (typename, dimensionality) for each variable of
	# operable types beginning with 'type'
	vars = dict()

	# Extract all of the variables of the managed types, so we can adjust
	# any required array sizes.
	for type in types:
		#re_var = re.compile ('\n\s*%%\(%s\)s.*::\s*([^()]*)([\s^\n]*\(([:,]*)\))*\n' % type)
		re_var = re.compile ('\n\s*%%\(%s\)s.*::\s*([^()]*)([\s^\n]*\(([:,]*)\))*.*\n' % type)

		for v in re_var.finditer(template):
			if v.group(2) == None:
				vars[v.group(1)] = (type, types[type])
			else:
				vars[v.group(1)] = (type, types[type] + v.group(3).count(":"))

			# If we are sizing an array based on the size of another, remove
			# this condition if the type is actually a scalar in this case.
			if types[type] == 0 and v.group(2) == None:
				re_fix = re.compile ('(\n\s*%%\(%s\)s.*::[\s^\n]*([^()\n]*))'
				                     '([\s^\n]*\(([^()\n]|\(([^()\n]|'
									 '\([^()\n]*\))*\))*\))*' % type)
				template = (template[0:v.start()] + 
				            re_fix.sub("\\1", template[v.start():], 1))

	# Replace the relevant occurrances
	for var in vars:
		# If we are considering dimensionality of 2 or larger
		if vars[var][1] > 1:
			dimstr = ":,"*(vars[var][1]-1)
			redim = '%s\s*\(' % var
			substr = ('%s(' % var) + dimstr
			re_redim = re.compile (redim)

			template = re_redim.sub (substr, template)

	return (template, config)

# Runtime entry point
if __name__ == '__main__':
	if len(sys.argv) < 2 or len(sys.argv) > 3:
		usage()
	else:
		fin = open(sys.argv[1], 'r')
		if len(sys.argv) == 2:
			fout = sys.stdout
		else:
			fout = open(sys.argv[2], 'w')

		print 'Input file: %s' % (fin.name)
		print 'Output file: %s' % (fout.name)

		config = read_config (fin)

		# Read in the template
		template = fin.read()
		(template,super_mod) = super_module(template, config)

		# Write the output file
		for s in config:
			# Adjust arrays to have the right number of dimensions
			(tmpl, cfg) = adj_arrays(template, config[s])

			fout.write(tmpl % cfg)
			fout.write("\n\n")
		fout.write(super_mod)

		fin.close()
		fout.close()
		

