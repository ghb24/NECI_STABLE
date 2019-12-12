#!/usr/bin/python

import sys

def extract_data(data_file):

    nsites = 0
    J = 0
    Ms = 0
    nbonds = 0
    bonds = {}

    f = open(data_file)

    first_line = True

    for line in f:
        values = line.split()
        if (not first_line):
            nbonds = nbonds + 1
            bonds[nbonds] = []
            bonds[nbonds].append(int(values[0]))
            bonds[nbonds].append(int(values[1]))
        if first_line:
            nsites = int(values[0])
            J = float(values[1])
            Ms = int(values[2]) 
            first_line = False

    f.close()

    return nsites, J, Ms, nbonds, bonds

def write_fcidump(nsites, J, Ms, nbonds, bonds):

    if (Ms + nsites) % 2 == 1:
        sys.exit("This Ms values is not possible for a lattice with this number of sites.")
    n_up = (Ms + nsites)/2

    print(' &FCI NORB=%d,NELEC= %d,MS2=  %d' % (2*nsites, 2*n_up, Ms))
    # Use sys.stdout.write instead of print, as we don't want a space after this.
    sys.stdout.write('  ORBSYM=')
    print '0,' * 2*nsites
    print('  ISYM=1 UHF=.TRUE.')
    print(' &END')
    for i in range(1, nbonds+1):
        print(' %20.13f%3d%3d%3d%3d' % (J/4, 2*bonds[i][0]-1, 2*bonds[i][0]-1, 2*bonds[i][1]-1, 2*bonds[i][1]-1))
        print(' %20.13f%3d%3d%3d%3d' % (-J/4, 2*bonds[i][0]-1, 2*bonds[i][0]-1, 2*bonds[i][1], 2*bonds[i][1]))
        print(' %20.13f%3d%3d%3d%3d' % (-J/4, 2*bonds[i][0], 2*bonds[i][0], 2*bonds[i][1]-1, 2*bonds[i][1]-1))
        print(' %20.13f%3d%3d%3d%3d' % (J/4, 2*bonds[i][0], 2*bonds[i][0], 2*bonds[i][1], 2*bonds[i][1]))
        print(' %20.13f%3d%3d%3d%3d' % (-J/2, 2*bonds[i][0]-1, 2*bonds[i][1]-1, 2*bonds[i][1], 2*bonds[i][0]))
        print(' %20.13f%3d%3d%3d%3d' % (-J/2, 2*bonds[i][0], 2*bonds[i][1], 2*bonds[i][1]-1, 2*bonds[i][0]-1))

if __name__ == '__main__':

    data_file = sys.argv[1]

    nsites, J, nbonds, Ms, bonds = extract_data(data_file)

    write_fcidump(nsites, J, nbonds, Ms, bonds)
