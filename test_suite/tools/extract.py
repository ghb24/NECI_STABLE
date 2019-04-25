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
# in the line where the data is stored. The fourth element, a logical,
# specifies whether the quantity should be repeated multiple times in the
# output file, if multiple simulations are being performed (i.e., replicas
# or excited states).
test_data = [
    ['energy_ground','GROUND', -1, False],
    ['ref_energy','<D0|H|D0>', -1, False],
    ['final_energy','Final energy estimate', -1, True],
    ['energy_rdm','*TOTAL ENERGY* CALCULATED USING THE *REDUCED DENSITY MATRICES*', -1, True],
    ['en2_energy','EN2 energy correction', -1, True],
    ['en2_tot_energy','*TOTAL ENERGY* including the EN2 correction', -1, True],
    ['max_error_rdm','MAX ABS ERROR IN HERMITICITY', 1, True],
    ['sum_error_rdm','SUM ABS ERROR IN HERMITICITY', 1, True],
    ['max_diff_trdm','MAX ABS DIFF IN HERMITICITY', 1, True],
    ['sum_diff_trdm','SUM ABS DIFF IN HERMITICITY', 1, True],
    ['1_rdm_diag_sum', 'SUM OF 1-RDM', -1, False],
    ['no_occ_sum', 'SUM OF THE N LARGEST NO OCCUPATION NUMBERS', -1, False],
    ['corr_entropy', 'CORRELATION ENTROPY  ', -1, False],
    ['hamil_sum','Sum of H', -1, False],
    ['overlap_sum','Sum of overlap matrix', -1, False],
    ['spec_low','Spectral weight at the lowest', -1, False],
    ['spec_high','Spectral weight at the highest', -1, False],
    ['ft_low','FT energy at lowest', -1, False],
    ['ft_high','FT energy at highest', -1, False],
    ['rep_est_var','Variational energy from replica_estimates', -1, True],
    ['rep_est_e_squ','Energy squared from replica_estimates', -1, True],
    ['rep_est_en2','EN2 estimate from replica_estimates', -1, True],
    ['rep_est_en2_new','New EN2 estimate from replica_estimates', -1, True],
    ['rep_est_precond','Preconditioned energy from replica_estimates', -1, True],
]

# The following are strings to be searched for which specify which simulation
# or estimate is about to be printed. For example, when using replica tricks
# or calculating excited states. The second index specifies where in the line
# the simulation/state label is specified.
simulation_labels = [
    ['Final energy estimate for state', 5],
    ['Stochastic error measures for RDM', -1],
    ['FINAL ESTIMATES FOR RDM', -1],
    ['PERFORMING ANALYSIS OF 1-RDM FOR STATE', -1],
    ['2-RDM ESTIMATES FOR STATE', -1],
    ['2-RDM ESTIMATES FOR TRANSITION', -4, -2, -1],
    ['REPLICA ESTIMATES FOR STATE', -1],
]

def extract_data(filename):
    '''Extract test data from the input file.'''

    names = []
    values = []
    sim_label_string = ''

    f = open(filename)

    for line in f:
        # Search for a string specifying which simulation/state the data will be for.
        for sim_string in simulation_labels:
            if sim_string[0] in line:
                words = line.split()
                if sim_string[0] == '2-RDM ESTIMATES FOR TRANSITION':
                    sim_label = words[sim_string[2]] + "_" + words[sim_string[1]] + words[sim_string[3]]
                    sim_label_string = "_" + sim_label
                else:
                    sim_label = words[sim_string[-1]]
                    sim_label = sim_label.rstrip(":")
                    sim_label = sim_label.rstrip("...")
                    sim_label_string = "_" + sim_label
        # Search for the data itself.
        for data in test_data:
            if data[1] in line:
                words = line.split()
                if data[3]:
                    # Use a header name with a simulation/state label.
                    names.append(data[0]+sim_label_string)
                else:
                    # Use a header name without a simulation/state label.
                    names.append(data[0])
                values.append(words[data[2]])

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
