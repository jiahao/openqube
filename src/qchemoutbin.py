#!/usr/bin/env python
"Read binary Q-Chem files"

import numpy
import os
FILE_SX      = 23
FILE_DENSITY_MATRIX=54
FILE_OVERLAP_MATRIX=320

def ReadOverlap(path = '', n = None):
    """
    Reads the overlap matrix.
    """
    filename = os.path.join(path, str(FILE_OVERLAP_MATRIX)+'.0')
    M = numpy.fromfile(filename, numpy.float64)
    if n is None:
        n = (M.size)**0.5
    M = M.reshape((n,n))
    return numpy.asmatrix(M)

def ReadDensity(path = '', n = None):
    """
    Reads the one-electron density matrix.

    Returns two numpy.matrices PA, PB
    for the alpha and beta density matrices respectively.
    """
    filename = os.path.join(path, str(FILE_DENSITY_MATRIX)+'.0')
    M = numpy.fromfile(filename, numpy.float64)
    if n is None:
        n = (M.size/2)**0.5
    PA, PB = M.reshape((2, n, n))
    return numpy.asmatrix(PA) + numpy.asmatrix(PB)

def ReadOverlapDeriv(path = ''):
    """
    Reads the overlap derivative matrix.
    """
    filename = os.path.join(path, str(FILE_SX)+'.0')
    M = numpy.fromfile(filename, numpy.float64)

    #TODO need to unpack sparse matrix format from Q-Chem.

    return M

if __name__ == '__main__':

    S = ReadOverlap('test')
    PA, PB = ReadDensity('test')
    print 'Overlap matrix'
    print S
    print 'Density matrix (alpha):'
    print PA
    print 'Density matrix (beta):'
    print PB
    print

    print 'Idempotency check: (should be 0.0)', numpy.linalg.norm(PA*S*PA - PA)

    print 'Alpha population', numpy.trace(S*PA)
    print 'Beta  population', numpy.trace(S*PB)

    Sx = ReadOverlapDeriv('test')
