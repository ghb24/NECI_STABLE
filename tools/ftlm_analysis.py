#!/usr/bin/python
'''ftlm_analysis.py [options] file

Calculate and output temperature-dependent energy results from a given set of
eigenvalues and corresponding initial overlaps. These should be output from
NECI from a KP-FCIQMC calculation.'''

import sys
import optparse
import pandas as pd
import math

def extract_data(data_files, cutoff):
    '''Extract the eigenvalues and initial overlaps from the data file
       provided.'''

    start_str = 'Eigenvalues and overlaps when keeping '+str(cutoff)

    f = open(data_file)
    have_data = False
    pairs = []
    for line in f:
        if not line.strip(): # If an empty line.
            have_data = False
        if have_data:
            values = line.split()
            pairs.append( [float(values[0]), float(values[1])] )
        if start_str in line:
            have_data = True
                
    return pairs

def calc_energy_contrib(pairs, minval, maxval, delta_beta):
    '''Calculate the contribution to the estimates of Tr(\rho*H) and Tr(\rho)
       from the eigenvalues and initial overlaps input in pairs, for the beta
       values calculated from minval, maxval and delta_beta.'''

    nbeta = int(math.ceil((maxval-minval)/delta_beta))+1

    beta_list = []
    energy_num_list = []
    trace_list = []

    for i in range(nbeta):
        beta = minval + i*delta_beta
        energy_num = 0.0
        trace = 0.0
        for [eigv, init_overlap] in pairs:
            energy_num += (init_overlap**2)*math.exp(-beta*eigv)*eigv
            trace += (init_overlap**2)*math.exp(-beta*eigv)
        beta_list.append(beta)
        energy_num_list.append(energy_num)
        trace_list.append(trace)

    results_dict = {'Beta' : pd.Series(beta_list), 
                    'Energy numerator' : pd.Series(energy_num_list),
                    'Trace' : pd.Series(trace_list)}

    results = pd.DataFrame(results_dict)

    return results

def parse_options(args):

    parser = optparse.OptionParser(usage = __doc__)
    parser.add_option('-m', '--min-plot', dest='minval', type='float', default=0.0,
                      help='The minimum beta to output results for.')
    parser.add_option('-n', '--max-plot', dest='maxval', type='float', default=10.0,
                      help='The maximum beta to output results for.')
    parser.add_option('-d', '--delta-beta', dest='delta', type='float', default=0.1,
                      help='The resolution in beta to plot.')
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

    # Loop over each data file (each of which contains data from one initial
    # configuration) and add in the contributions from each of those files.
    for data_file in data_files:
        # Extract the eigenvalues and overlaps with the intial Krylov vector.
        pairs = extract_data(data_file, options.cutoff)
        results = calc_energy_contrib(pairs, options.minval, options.maxval, options.delta)

        try:
            final_results[['Energy numerator', 'Trace']] = \
                final_results[['Energy numerator', 'Trace']] + results[['Energy numerator', 'Trace']]
        except NameError:
            # This is entered the first time, when final_results doesn't yet
            # exist.
            final_results = results

    print final_results.to_string(index=False)
