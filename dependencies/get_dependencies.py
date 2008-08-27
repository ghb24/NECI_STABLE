#!/usr/bin/python

import re,sys

utilities='bookkeeper_pointers common_routines default_sets input Logging memory_manager MemoryManager parameters precision record_handler record_handler_arrays record_handler_pointers run_data'.split()

is_utility=re.compile('\\b%s\\b' % ('\\b|\\b'.join(utilities)),re.I)
module_regex=re.compile('(^ *module )(\\b[a-z_]+\\b)',re.I)
end_module_regex=re.compile('end * module',re.I)
use_regex=re.compile('(^ *use )(\\b[a-z_]+\\b)',re.I)

class use_sets(object):
    def __init__(self):
        self.main=set()
        self.utilities=set()

def update_use_list(using_object,used_object,use_list):
    dependency=re.sub('\.','_','%s -> %s' % (using_object,used_object))
    if re.match(is_utility,using_object) or re.match(is_utility,used_object):
        use_list.utilities.update([dependency])
    else:
        use_list.main.update([dependency])

def output_dependencies(dep_name,dependencies):
    f=open('%s.dot' % dep_name,'w')
    f.write('digraph %s {\n' % dep_name)
    for line in dependencies:
        f.write('\t%s;\n' % line)
    f.write('}')
    f.close()

modules_dependencies=use_sets()
files_dependencies=use_sets()

file_list=sys.argv[1:]

for file in file_list:
    is_module=False
    source_code=open(file,'r').readlines()
    for line in source_code:
        line_search=re.findall(module_regex,line)
        if line_search:
            module_name=line_search[0][1]
            is_module=True
        elif re.search(end_module_regex,line):
            is_module=False
        else:
            line_search=re.findall(use_regex,line)
            if line_search:
                used_object=line_search[0][1]
                if is_module:
                    update_use_list(module_name,used_object,modules_dependencies)
                else:
                    update_use_list(file,used_object,files_dependencies)

output_dependencies('modules_use_main',modules_dependencies.main)
output_dependencies('modules_use_utilities',modules_dependencies.utilities)
output_dependencies('files_use_main',files_dependencies.main)
output_dependencies('files_use_utilities',files_dependencies.utilities)
