#!/usr/bin/env python3

'''extract_eigv.py -l nlowdin file.*

Extract the Hamiltonian eigenvalues from file.* files and move them to a
file with the same form as EIGV_DATA files output by NECI. The 'file.*' files
should be called file.I, where I is an iteration label which will be used in
the EIGV_DATA file, and should have the format of lowdin.* files output by
NECI when running KP-FCIQMC calculations.'''

import optparse
import os
import sys

def extract_data(data_file, target_nlowdin):
    '''Extract the eigenvalues from data_file, using a Lowdin cutoff of
       nlowdin_target. If this Lowdin cutoff does not exist then use the
       largest cutoff that does.'''

    overlap_start_str = 'Overlap matrix eigenvalues'
    start_str = 'Eigenvalues and overlaps when keeping ' + str(target_nlowdin)

    nlowdin = 0
    have_overlap = False
    have_data = False
    energies = []

    f = open(data_file)

    for line in f:
        if not line.strip(): # If an empty line.
            # If we've just finished looking at overlap eigenvalues.
            if have_overlap:
                start_str = 'Eigenvalues and overlaps when keeping ' + str(min(nlowdin, target_nlowdin))

            have_overlap = False
            have_data = False

        if have_overlap:
            # Count the number of positive overlap eigenvalues.
            values = line.split()
            overlap_eval = float(values[0])
            if overlap_eval > 0.0:
                nlowdin += 1

        if have_data:
            values = line.split()
            energies.append(float(values[0]))

        if overlap_start_str in line:
            have_overlap = True
        if start_str in line:
            have_data = True

    return nlowdin, energies

def parse_options(args):
    '''Parse arguments and files from the command line.'''

    parser = optparse.OptionParser(usage = __doc__)
    parser.add_option('-l', '--lowdin-cutoff', dest='cutoff', type='int', default=5,
                      help='The number of eigenvectors to keep in the Lowdin '
                      'orthogonalisation procedure.')
    (options, filenames) = parser.parse_args(args)

    if len(filenames) == 0:
        parser.print_help()
        sys.exit(1)

    return (options, filenames)

if __name__ == '__main__':

    (options, data_files) = parse_options(sys.argv[1:])

    # Get the stem and extension (iteration label, together with a '.').
    stem, iteration = os.path.splitext(data_files[0])

    # Find all the iterations labels and then sort them.
    iterations = []
    for data_file in data_files:
        iterations.append(int(os.path.splitext(data_file)[1][1:]))
    iterations.sort()

    # Print the header.
    sys.stdout.write('# 1. Iteration')
    icolumn = 1
    for ivec in range(1,options.cutoff+1):
        icolumn += 1
        header = str(icolumn) + '. Energy ' + str(ivec)
        sys.stdout.write('%22s' % header)
    sys.stdout.write('\n')

    # Extract and print information for each iteration.
    for iteration in iterations:
        filename = stem + '.' + str(iteration)
        nlowdin, energies = extract_data(filename, options.cutoff)
        sys.stdout.write('     %9d' % iteration)
        for energy in energies:
            sys.stdout.write('   %19.12e' % energy)
        if nlowdin < options.cutoff:
            for i in range(nlowdin+1, options.cutoff+1):
                sys.stdout.write('            NaN       ')
        sys.stdout.write('\n')
