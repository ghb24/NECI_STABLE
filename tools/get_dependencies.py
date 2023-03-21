#!/usr/bin/env python
'''
Analyse Fortran source code and output (in the dot language) the dependencies on modules in the source code.

Usage:

get_dependencies.py [files containing source code]

Output in the dot language:
* modules depending on other modules.
* non-modules files depending on modules.
These are further split up into dependencies involving utilty or data modules (which must be specified inside get_dependencies).

The resultant *.dot files are in the output_dir directory, as specified in the
configuration section.  Graphs can be produced using a suitable parser.

For more information on the dot langauge, see http://www.graphviz.org.
'''

from __future__ import print_function
import re,sys

__author__='James Spencer'

#----------------------------------------------------------------------

# Configuration.

# List of utilities modules and data modules (splitting a string is more
# convenient than actually specifying a list of strings).  
# Determine whether a module is a utility or data module via regular
# expressions.  Case is ignored.
# A dependency (in either direction) that involves one of the modules
# listed is placed in the graph of that category.
# All other modules are put into the "main" category.
utilities_modules='bookkeeper_pointers common_routines default_sets global_utilities helem input Logging memory_manager MemoryManager mpi parameters precision record_handler record_handler_arrays record_handler_pointers run_data timing'.split()
data_modules='CalcData CPMDData IntegralsData SymData SystemData'.split() 

# Output directory for *.dot files containing dependencies. (Must exist.)
output_dir='dependencies'

#----------------------------------------------------------------------

# Regular expressions for module types.
# Warning: if used as a module and the user changes the modules lists, these
# are *not* updated.
is_utility_object=re.compile('\\b%s\\b' % ('\\b|\\b'.join(utilities_modules)),re.I)
is_data_object=re.compile('\\b%s\\b' % ('\\b|\\b'.join(data_modules)),re.I)

# Regular expressions for the start, end and usage of modules.
# Note that end of module detection relies on using 'end module' rather
# than just 'end'.
module_regex=re.compile('(^ *module )(\\b[a-z_]+\\b)',re.I)
end_module_regex=re.compile('end * module',re.I)
use_regex=re.compile('(^ *use )(\\b[a-z_]+\\b)',re.I)

class use_sets(object):
    '''Store the dependency entries for dependencies on data, utility and all other modules in one handy object.'''
    def __init__(self):
        self.data=set()
        self.main=set()
        self.utilities=set()
    def __repr__(self):
        s='data modules dependencies: %s\nutilties modules dependencies: %s\nmain modules dependencies: %s' % (repr(self.data),repr(self.utilities),repr(self.main))
        return s

def update_use_list(using_object,used_object,use_list):
    '''Add dependency (in dot format) to the appropriate element of the use_list (itself a use_sets object).'''
    dependency=re.sub('\.','_','%s -> %s' % (using_object,used_object))
    if re.match(is_data_object,using_object) or re.match(is_data_object,used_object):
        use_list.data.update([dependency])
    elif re.match(is_utility_object,using_object) or re.match(is_utility_object,used_object):
        use_list.utilities.update([dependency])
    else:
        use_list.main.update([dependency])

def output_dependencies(dep_name,dependencies):
    '''Write out dependencies to files.'''
    f=open('%s/%s.dot' % (output_dir,dep_name),'w')
    f.write('digraph %s {\n' % dep_name)
    for line in dependencies:
        f.write('\t%s;\n' % line)
    f.write('}')
    f.close()

def main(file_list):
    '''Generate and output dependencies of the Fortran source code files in the file list.'''
    modules_dependencies=use_sets()
    files_dependencies=use_sets()

    for file in file_list:
        is_module=False
        source_code=open(file,'r').readlines()
        for line in source_code:
            line_search=re.findall(module_regex,line)
            if line_search:
                # We convert all module names to Title_Case to avoid 
                # any issues with different cases being used.
                module_name=line_search[0][1].title()
                is_module=True
            elif re.search(end_module_regex,line):
                is_module=False
            else:
                line_search=re.findall(use_regex,line)
                if line_search:
                    used_object=line_search[0][1].title()
                    if is_module:
                        update_use_list(module_name,used_object,modules_dependencies)
                    else:
                        update_use_list(file,used_object,files_dependencies)

    output_dependencies('modules_use_data',modules_dependencies.data)
    output_dependencies('modules_use_main',modules_dependencies.main)
    output_dependencies('modules_use_utilities',modules_dependencies.utilities)
    output_dependencies('files_use_data',files_dependencies.data)
    output_dependencies('files_use_main',files_dependencies.main)
    output_dependencies('files_use_utilities',files_dependencies.utilities)

if __name__=='__main__':
    if len(sys.argv)==1 or sys.argv[1]=='-h' or sys.argv[1]=='--h' or sys.argv[1]=='-help' or sys.argv[1]=='--help':
        print(__doc__) # Wow.  Magic...
        sys.exit()
    else:
        file_list=sys.argv[1:]
        main(file_list)
