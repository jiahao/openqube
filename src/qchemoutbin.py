#!/usr/bin/env python
"Read binary Q-Chem files"

import numpy
import os

class QChemBinaryOutput:
    "Read binary Q-Chem files"

    #Enumerated file numbers in Q-Chem
    FILE_SX      = 23
    FILE_DENSITY_MATRIX=54
    FILE_OVERLAP_MATRIX=320

    def __init__(self, path = ''):
        self.path = path

    def OverlapMatrix(self, n = None):
        """
        Reads the overlap matrix.
        """
        filename = os.path.join(self.path, str(self.FILE_OVERLAP_MATRIX)+'.0')
        M = numpy.fromfile(filename, numpy.float64)
        if n is None:
            n = (M.size)**0.5
        M = M.reshape((n,n))
        return numpy.asmatrix(M)
    
    def DensityMatrix(self, SeparateSpin = False, n = None):
        """
        Reads the one-electron density matrix.
    
        Returns two numpy.matrices PA, PB
        for the alpha and beta density matrices respectively.
        """
        filename = os.path.join(self.path, str(self.FILE_DENSITY_MATRIX)+'.0')
        M = numpy.fromfile(filename, numpy.float64)
        if n is None:
            n = (M.size/2)**0.5
        DensityMat = M.reshape((2, n, n))
        if SeparateSpin:
            return numpy.asmatrix(DensityMat[0,:,:]), numpy.asmatrix(DensityMat[1,:,:])
        else:
            return numpy.asmatrix(DensityMat[0,:,:] + DensityMat[1,:,:])

    def OverlapMatrixDerivative(self):
        """
        Reads the overlap derivative matrix.
        @todo need to unpack sparse matrix format from Q-Chem.
        """
        filename = os.path.join(self.path, str(self.FILE_SX)+'.0')
        M = numpy.fromfile(filename, numpy.float64)
    
        return M

if __name__ == '__main__':

    QBin = QChemBinaryOutput('test')
    S = QBin.OverlapMatrix()
    PA, PB = QBin.DensityMatrix(SeparateSpin = True)
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

    Sx = QBin.OverlapMatrixDerivative()

