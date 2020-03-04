import os

def writeDMRGconfig(csf, savePath=os.getcwd()):
    '''for easy comparison of my results with Sandeeps Block code execute CSHOF for given CSF + connected CSFs'''
    assert len(csf) > 0 and len([n for n in csf if 0 <= n <= 3]) == len(csf),'no proper csf given!'

    nOrbitals, nElectrons, spin = calcElectronConfiguration(csf);

    dmrgConfig = open(savePath + os.path.sep + 'dmrg.conf', 'w');

    # write configuration specific infos
    dmrgConfig.write('nelec ' + str(nElectrons) + '\n');
    dmrgConfig.write('spin ' + str(spin) + '\n');
    
    # write standard rest:
    dmrgConfig.write('\nschedule\n0 50 1.0e-5 10.0\n1 50 1.0e-6 1.0e-4\n3 50 1.0e-6 1.0e-5\n4 50 1e-6 0.0\nend\n');
    dmrgConfig.write('\nmaxiter 7\ntwodot_to_onedot 5\nsweep_tol 1e-10\n\norbitals FCIDUMP\nsym d2h\n\n\n');
    dmrgConfig.write('warmup wilson\n\nhf_occ integral\n\nprefix scratch/\nscreen_tol 0.0\n\n\nnoreorder\noutputlevel 5\n');

    dmrgConfig.close();

def writeFCIDUMP(csf, i, j, k = 0, l = 0, savePath=os.getcwd()):
    '''write FCIDUMP file for test of single and double excitations with sandeeps code'''
    assert len(csf) > 0 and len([n for n in csf if 0 <= n <= 3]) == len(csf),'no proper csf given!'
    assert isinteger(i) and isinteger(j) and isinteger(k) and isinteger(l)
    assert i > 0 and j > 0 and k >= 0 and l >= 0

    nOrbitals, nElectrons, spin = calcElectronConfiguration(csf);

    FCIDUMP = open(savePath + os.path.sep + 'FCIDUMP','w');
    
    # header part:
    FCIDUMP.write('&FCI NORB= ' + str(nOrbitals) + ',NELEC= ' + str(nElectrons) + ',MS2= ' + str(spin) + ',\n');
    FCIDUMP.write('ORBSYM=' + nOrbitals*'1,' + '\nISYM=1,\n&END\n');

    # overlap integral:
    FCIDUMP.write('1.0   ' + str(i) + '   ' + str(j) + '   ' + str(k) + '   ' + str(l) + '\n');
    FCIDUMP.write('0.0   0   0   0   0');

    FCIDUMP.close();

  
def writeDeterminants(csf, connectedCSFs, savePath=os.getcwd()):
    '''write determinants input file for Block code for given connected CSF to compare results with sandeeps code'''
    
    # pretend parent csf to connectedCSF list:
    connectedCSFs.reverse();
    connectedCSFs.append(csf);
    connectedCSFs.reverse();

    bitString = convertStepVectorToBitString(connectedCSFs);

    determinants = open(savePath + os.path.sep + 'determinants','w');

    determinants.write(str(len(bitString)) + '\n' );
    
    for strings in bitString:
        determinants.write(strings + '\n');

    determinants.close();


def convertStepVectorToBitString(csfs):
    '''converts a given collection of csfs in stepVector representation to string in binary format to use for sandeeps
    Block code'''

    out = [];

    for item in csfs:
        string = '';
        for step in item:
            if step == 0:
                string += '0 0  ';

            elif step == 1:
                string += '1 0  ';

            elif step == 2:
                string += '0 1  ';

            elif step == 3:
                string += '1 1  ';

        out.append(string);

    return out

                
                

def calcElectronConfiguration(csf):
    ''' calc. the electron configuration for a given csf'''
    assert len(csf) > 0 and len([n for n in csf if 0 <= n <= 3]) == len(csf),'no proper csf given!'

    nOrbitals = len(csf);
    nElectrons = sum([2 for n in csf if n == 3]) + sum([1 for n in csf if 1 <= n <= 2]);
    spin = sum([1 for n in csf if n == 1]) + sum([-1 for n in csf if n == 2]);

    return nOrbitals, nElectrons, spin



