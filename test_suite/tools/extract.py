#!/usr/bin/env python
'''extract.py [options] file

Extract test suite data and output it in a format which can be used by
testcode2.'''

import optparse
import sys

# Test data to be searched for.
# Each element in this list is a list with information about a single test
# value. The first element gives the name, the second a string that should be
# searched for in the file, and the third gives the position of the data value
# in the line where the data is stored.
test_data = [
    ['energy_ground','GROUND', -1],
    ['ref_energy','<D0|H|D0>', -1],
    ['energy_summed','Summed', -1],
    ['energy_rdm','*TOTAL ENERGY* CALCULATED USING THE REDUCED DENSITY MATRICES', -1],
    ['max_error_rdm','MAX ABS ERROR IN HERMITICITY', 1],
    ['sum_error_rdm','SUM ABS ERROR IN HERMITICITY', 1],
    ['overlap_sum','Sum of H', -1],
    ['hamil_sum','Sum of overlap matrix', -1],
    ['spec_low','Spectral weight at the lowest', -1],
    ['spec_high','Spectral weight at the highest', -1],
    ['ft_low','FT energy at lowest', -1],
    ['ft_high','FT energy at highest', -1]
]

def extract_data(filename):
    '''Extract test data from the input file.'''

    names = []
    values = []

    f = open(filename)

    for line in f:
        for data in test_data:
            if data[1] in line:
                words = line.split()
                names.append(data[0])
                values.append(words[data[-1]])

    f.close()

    return names, values

def write_data(names, values, padding):
    '''Write test data in a clean format.'''

    # Print data names.
    for (name, value) in zip(names, values):
        format = '%%-%is' % (max(len(value), len(name))+padding)
        sys.stdout.write(format % name)

    sys.stdout.write('\n')

    # Print the data values themselves.
    for (name, value) in zip(names, values):
        format = '%%-%is' % (max(len(value), len(name))+padding)
        sys.stdout.write(format % value)

    sys.stdout.write('\n')

def parse_options(args):
    '''Read and return filename and any options present.'''

    parser = optparse.OptionParser(usage = __doc__)
    (options, filename) = parser.parse_args(args)

    if len(filename) != 1:
        parser.print_help()
        sys.exit(1)

    return filename[0]

if __name__ == '__main__':
    filename = parse_options(sys.argv[1:])

    names, values = extract_data(filename)
    write_data(names, values, padding=2)
