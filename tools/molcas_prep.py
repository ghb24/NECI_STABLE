#!/usr/bin/python -tt
"""
Produce a directory of output files suitable for incorporation into molcas
from the NECI source code

Usage:
    molcas_prep.py [NECI directory] [tgt_directory]

Note - this script creates the target directory, which must not exist prior to
running it. This protects against accidental overwriting of files.

--> Runs the templating engine to produce .F90 files from the .F90.template files
--> Renames all .F* files to .f* for the molcas build system
--> Ensures that each file contains only one module. This requires splitting files
    up.
"""

import os
import sys
import shutil
import re
import f90_template


def usage():
    """
    Prints the usage description for the program
    """
    print __doc__


def file_error():
    """
    Give a nice looking error if something happens...
    """
    print "A file-handling error occurred"
    sys.exit(-1)


def process_f(dir, fn, tgt_dir, tmp_dir):
    """
    Process a .F file for molcas-suitability
    """
    # Get the source and target full filenames. Note we convert .F to .f
    root, ext = os.path.splitext(fn)
    src_file = os.path.join(dir, fn)
    tgt_file = os.path.join(tgt_dir, "{}.f".format(root.lower()))
    shutil.copyfile(src_file, tgt_file)


def process_f90(dir, fn, tgt_dir, tmp_dir):
    """
    Process an F90 file for molcas-suitability. Note that F90 files may contain modules, which
    makes them different to F files above.
    """
    # Get the source and target full filenames. Note we convert .F90 to .f90
    root, ext = os.path.splitext(fn)
    src_file = os.path.join(dir, fn)
    tgt_file = os.path.join(tgt_dir, "{}.f90".format(root.lower()))
    assert ext == ".F90"
    assert src_file != tgt_file

    # We wish to test how many modules are in the specified file.
    re_mod = re.compile('(\s*module(.|\n)*?\n\s*end module)', re.IGNORECASE)
    with open(src_file, 'r') as fin:
        contents = fin.read()
        found_mods = re_mod.findall(contents)

        # If there is only one module (or fewer) in the file, then we can just copy this
        # file directly.
        if len(found_mods) in (0, 1):
            shutil.copyfile(src_file, tgt_file)
        else:
            # We wish to move all modules other than the _last_ module into other
            # files.
            print "Rejecting F90 file {} with {} modules".format(src_file, len(found_mods))
            return


def process_f90_template(dir, fn, tgt_dir, tmp_dir):
    """
    Process an F90 template file
    """
    # Get the equivalent F90 filename
    root, ext = os.path.splitext(fn)
    src_file = os.path.join(dir, fn)
    tmp_file = os.path.join(tmp_dir, root)

    # Open these files for use, and call the template processing utility.
    # We can use this directly without calling f90_template.py as a script.
    with open(src_file, 'r') as fin:
        with open(tmp_file, 'w') as fout:
            filelist = f90_template.process_file(fin, fout, silent=True, multifile=True)

    # Now we process this as though it were a normal F90 file
    for written_fn in filelist:
        name = os.path.basename(written_fn)
        process_f90(tmp_dir, name, tgt_dir, tmp_dir)


def file_direct_copy(dir, fn, tgt_dir, tmp_dir):
    """
    Directly copy the specified file into the target directory
    """
    shutil.copyfile(os.path.join(dir, fn), os.path.join(tgt_dir, fn.lower()))


def drop_file(dir, fn, tgt_dir, tmp_dir):
    """
    The specified file should not be used!
    """
    pass


def process_files(src_dir, tgt_dir, tmp_dir):
    """
    Do the actual file processing
    """

    # This is the mapping for how the files ought to be treated
    file_map = {
        'F': process_f,
        'F90': process_f90,
        'F90.template': process_f90_template,
        'cpp': file_direct_copy,
        'h': file_direct_copy
    }

    # Ensure that the target and tmp directories exist
    os.mkdir(tgt_dir)
    os.mkdir(tmp_dir)

    # Ensure that the target directory is clean

    # List the files present in the directory (be recursive)
    for root, dirs, files in os.walk(src_dir, onerror=file_error):
        for fn in files:
            try:
                # Determine what to do by the files extension. If unknown, then
                # discard the file.
                ext = fn.split('.', 1)[1]
                file_map.get(ext, drop_file)(root, fn, tgt_dir, tmp_dir)
            except IndexError:
                # Files with no index should be dropped
                drop_file(root, fn, tgt_dir, tmp_dir)

    # Remove the temporary files directory
    shutil.rmtree(tmp_dir)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        usage()
        sys.exit(-1)

    src_dir = os.path.join(sys.argv[1], 'src')
    tgt_dir = sys.argv[2]
    tmp_dir = os.path.join(tgt_dir, 'tmpfiles')
    print "Source directory: ", src_dir
    print "Target directory: ", tgt_dir
    print "Temporary directory: ", tmp_dir

    # And kick off the calculation
    process_files(src_dir, tgt_dir, tmp_dir)


