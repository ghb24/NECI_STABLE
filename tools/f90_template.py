#!/usr/bin/python -tt
'''Produce a .f90 file from the given .F90.template file to simulate the
effects of using templates/generic functions in a language such as c++

Usage:
	f90_template.py [infile] [outfile]

The template files are split into two sections. The first is a config section.
Each desired configuration is demarcated with a section header of the type:

[section_name]

Beneath this are name/value pairs. These are used as a literal substitutions
into the second section. This section implements a form of inheritance. All
values declared in a section are available to the following sections unless
explicitly changed.

Note that beginning a type-name label with "type" implements a couple of
additional features.

************************

The second section implements a module. Currently the script only supports 
using one module. Any of the labels declared in the first section are
available, and may be used as follows:

	%(label-name)s

In addition, there are a couple of extra labels present:

	%(dim-name)s - For an array called 'name', this contains the dimenionality
	               of the array.
				   This is a little bit fragile at the moment, so should be 
				   used with a bit of caution.
	%(name)s - Contains the current configuration name as given in the first
	           section.

If a type has been declared in the first section, and that type contains a
statement of dimension(...), then the script is able to adjust the
dimensionality of of any references to that array so that array slices of the
correct size are used. Please see src/lib/quicksort.F90.template for an
example.

An functions must be declared using the result(name) structure, i.e.

	function test_fn (args) result(res)

		type :: res
		...
	end function

These functions, as well as all subroutines will be renamed to avoid name
space clashes between templated modules, and an interface generated to allow access as expected.

************************

If the property 'conditional_enable' is defined, it is only valid for the
particular configuration it is specified in.

It will be the condition placed in:

	#if *
	#endif

block around the configuration.

************************

Supermodule:
Once multiple, subtly renamed, modules have been produced in the output file,
they are then all included into an overall module, with the specified name, 
which can be used in the code.

If there is anything to include in the module which does not need to be
templated, then it can be included here. See allocate_shared.F90.template.
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
		if line and line[0] == '=':
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
		config[s] = dict()
		if sprev is not None:
			config[s].update(config[sprev])
		config[s].update({'conditional_enable' : None})
		config[s].update(parser.items(s))
		config[s].update({'name' : s})

		if sprev:
			sys.stdout.write(", ")
		sys.stdout.write("%s" % (s))
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
	if not m:
	   print("Did not find module string. Exiting")
	   exit() 
	print 'Generating super module: %s' % m.group(2)

	# Change module in template to be a sub-module with associated name
	substr = "%s_%%(name)s\n" % m.group(1)
	template = re_mod.sub (substr, template)

	# Extract any lines for the super module from the template file
	re_supermod = re.compile ('(\n\s*supermodule[\s^\n]*([^\s\n]*))\n((.|\n)*?)(end supermodule)')
	m_super = re_supermod.search(template)
	if m_super:
		print 'Found super module: %s.' % m.group(2)
		template = template[:m_super.start()] + template[m_super.end()+1:]

	# Construct super module
	super_mod = m.group(1) + "\n"
	for s in config:
		if config[s]['conditional_enable'] is not None:
			super_mod += '#if %s\n' % config[s]['conditional_enable']
		super_mod += "    use " + m.group(2) + "_" + s + "\n"
		if config[s]['conditional_enable'] is not None:
			super_mod += '#endif\n'

	if m_super: super_mod += m_super.group(3)
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
	
	# Vars contains tuples of (typename, dimensionality, dimstring) for each variable of
	# operable types beginning with 'type'.  dimstring is the string that was used for declaring the variable's dimensions

	# We split the template into sections for each function or subroutine
	tsplitter=re.compile("\n[^!\n]*(subroutine|function)",flags=re.IGNORECASE)

	# Locate any local interfaces.
	re_interfaces = re.compile("(\n[^!\n]*interface)((.|\n)+?)(\n[^!\n]*end\s+interface)")

	# Find splitting points for subroutines excluding those within
	# interface statements.
	split_locs = []
	for s in tsplitter.finditer(template):
		bInterfaced = False
		for i in re_interfaces.finditer(template):
			if s.start() > i.start() and s.end() < i.end():
				bInterfaced = True
				break
		if not bInterfaced:
			split_locs.append(s.end()-2)

	newtemplate=template[:split_locs[0]]
	for i in range(len(split_locs)):
		templpart = ""
		if (i == len(split_locs)-1):
			templpart = template[split_locs[i]:]
		else:
			templpart = template[split_locs[i]:split_locs[i+1]]

		vars = dict()

		# Extract all of the variables of the managed types, so we can adjust
		# any required array sizes.

		# types is a dict with keys being the type line from config 
		# (e.g. type1 = integer, dimension(:)) and the values being the
		# number of dimensions
		for type in types:
				#re_var = re.compile ('\n\s*%%\(%s\)s.*::\s*([^()]*)([\s^\n]*\(([:,]*)\))*\n' % type)
				re_var = re.compile ('\n\s*%%\(%s\)s.*::\s*([^()]*)([\s^\n]*\(([:,]*)\))*.*\n' % type)
				
				# Search for lines of the form
				# TYPE :: VARIABLE
				# TYPE :: VARIABLE(:)
				# TYPE :: VARIABLE(:,:) etc.
				
				v = re_var.search(templpart)
				offset = 0
				while v:
					if v.group(2) == None:
						vars[v.group(1)] = (type, types[type],"")
					else:
						vars[v.group(1)] = (type, types[type] + v.group(3).count(":"),v.group(3))

					# If we are sizing an array based on the size of another, remove
					# this condition if the type is actually a scalar in this case.
					if types[type] == 0 and v.group(2) == None:
						re_fix = re.compile ('(\n\s*%%\(%s\)s.*::[\s^\n]*([^\(\)\n]*))'
											 '([\s^\n]*\(([^\(\)\n]|\(([^\(\)\n]|'
											 '\([^\(\)\n]*\))*\))*\))*' % type)
						templpart = (templpart[0:v.start()+offset] + 
						re_fix.sub("\\1", templpart[v.start()+offset:], 1))
					offset = offset + v.start() + 1
					v = re_var.search(templpart[offset:])
		
		# vars is a dict with keys being the name of the variable and values 
		# being (TYPE,TOTALDIMS,DIMSTRING)
		# TOTALDIMS is the sum of dims in the TYPE and of the VARIABLE
		# DIMSTRING is the string (e.g. :,: ) which was used to declare the 
		# variable's dimensions

		# Replace the relevant occurrances
		for var in vars:
			config['dim-%s' % var] = '%d' % vars[var][1]
			# Only consider types with a dimenison() term to be adjustable
			# based on the num of dimensions in type. Need to get this far
			# to set dim-%s.
			if types[vars[var][0]] != 0:  # If the type is not a scalar
				# If we are considering dimensionality of 2 or larger
				dimstr = ":,"*(vars[var][1]-1)
				#redim = '([\( ,]+)%s\s*\(' % (var)
				#substr = ('\\1%s(' % var) + dimstr
				redim = '([\( ,=])%s\s*\(' % (var)
				substr = ('\\1%s(' % var) + dimstr

				# If our original dimension string was empty (i.e. didn't 
				# contain anything special we had to leave) then we need to
				# include the full number of dimensions, not one fewer
				#  - we do this by changing VARIABLE() to VARIABLE(:)

				if vars[var][2]=="":
					re_redim = re.compile (redim+"\)")
					templpart = re_redim.sub (('\\1%s(' % var)+":)", templpart)
				re_redim = re.compile (redim)
				templpart = re_redim.sub (substr, templpart)

			# If we've ended up with a scalar, just remove the (:)
			elif vars[var][1]==0:
				redim = '%s\s*\(%s\)' % (var,vars[var][2])
				re_redim = re.compile(redim)
				templpart = re_redim.sub(var,templpart)
				
		newtemplate += templpart 


	# Now we want to search through for any interfaces. We then want to remove the (*)
	# from any arrys which have them.
	split_ints = []
	last_end = 0
	new_template2 = ""
	for s in re_interfaces.finditer (newtemplate):

		# Where does this interface start/end
		start, end = s.start(), s.end()
		templpart = template[start:end+1]

		# Fill in any unchanged region.
		new_template2 += template[last_end:start]

		# Find and fix any variables of the managed types.
		for type in types:

			# If the type is a scalar, then we need to remove the assumed-type part
			# of it
			if types[type] == 0:

				re_assumed = re.compile ('(\n\s*%%\(%s\)s.*::\s*[^()]*)([\s^\n]*\(\*\))(.*\n)' % type)

				v = re_assumed.search(templpart)
				offset = 0
				while v:

					# Fix this section by removing assumed-type parts
					templpart = (templpart[0:v.start()+offset] +
						re_assumed.sub("\\1\\3", templpart[v.start()+offset:], 1))

					offset = offset + v.start() + 1
					v = re_assumed.search(templpart[offset:])


		# We want to start after this region.
		new_template2 += templpart
		last_end = end

	# Fill in the rest
	new_template2 += newtemplate[last_end:]

	if config['conditional_enable'] is not None:

		new_template2 = '#if %s\n' % config['conditional_enable'] + new_template2 + '\n#endif\n'

	return (new_template2, config)




def interface_procs (template):
	'''Create module procedure interfaces for all of the module procedures'''

	print "Interface generation"

	# Do we have procedures in this template file?
	re_contains = re.compile ('\n\s*contains\s*\n',flags=re.IGNORECASE)
	re_end_mod = re.compile ('\n\s*end\s*module')
	m = re_contains.search(template)
	if m:
		interface = "\n"

		# Find all of the functions or subroutines. Note that it requires functions to
		# be new style functions, with results declared in result(ret) format.
		#re_proc = re.compile ('\n\s*(subroutine|function)[\s^\n]+\(\w+\)[\s^\n]+\(')
		re_proc = re.compile ('(\n\s*((elemental|pure)[\s^\n]+)*(subroutine|function)[\s^\n]+)(\w+)([\s^\n]*\()',re.I)

		offset = m.start()
		m_end_mod = re_end_mod.search(template[offset:])
		proc = re_proc.search(template[offset:])
		while proc:
			print "Procedure: ",proc.group(5)
			# For each procedure, append _%(name)s to all of the names.
			template = (template[0:offset + proc.start()] +
					   re_proc.sub("\\1\\5_%(name)s\\6", template[proc.start()+offset:], 1))
			if proc.start() < m_end_mod.end():
				interface += ('    interface %s\n'
							  '        module procedure %s_%%(name)s\n'
							  '    end interface\n' % (proc.group(5), proc.group(5)))
			offset = offset + proc.end() + 10
			m_end_mod = re_end_mod.search(template[offset:])
			proc = re_proc.search(template[offset:])

		interface += "contains\n"
		template = template[0:m.start()] + interface + template[m.end():]

	return template

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
		template = interface_procs (template)

		# Write the output file
		for s in config:
			# Adjust arrays to have the right number of dimensions
			(tmpl, cfg) = adj_arrays(template, config[s])
			try:
				fout.write(tmpl % cfg)
			except ValueError, e:
				print 'Value error: ', e
				fout.write(tmpl)
			except Exception, e:
				print tmpl
				print "Exception: ", e
				raise
			fout.write("\n\n")
		fout.write(super_mod)

		fin.close()
		fout.close()
		

