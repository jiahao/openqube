#!/usr/bin/env python
import os
from numpy import real
from numpy.linalg import eigh
def WriteDensityToMolden(filename, P):
    """This function overwrites an existing Molden file with a new density
    matrix P. This is accomplished by diagonalizing the density matrix P
    and writing its natural orbitals to the file"""

    assert os.path.exists(filename), 'Need existing Molden file as template'

    #Calculate (nonorthogonal) natural orbitals and their occupation numbers
    [occnos, NMOs] = eigh(P)

    buf = []
    orbidx = P.shape[0] #Orbital index
    state = 'seek'
    for line in open(filename):
        newline = line[:-1]
        if state == 'seek':
            if '[MO]' in line:
                state = 'MO'
        elif state == 'MO':
            if 'Ene=' in line:
                #XXX What are the energies of natural orbitals, really?
                orbidx -= 1
            elif 'Occup=' in line:
                newline = 'Occup=%f' % occnos[orbidx]
            else:
                try:
                    t = line.split()
                    aoidx = int(t[0])-1
                    float(t[1])
                    newline = t[0] + ' ' + str(real(NMOs[aoidx, orbidx]))
                except (ValueError, IndexError):
                    pass

        else:
            raise ValueError, 'Unknown state'
        buf.append(newline)

    #Now write modified Molden file
    f = open(filename, 'w')
    f.write('\n'.join(buf))
    f.close()
