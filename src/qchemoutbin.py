#!/usr/bin/env python
"Read binary Q-Chem files"

import numpy, os

class QChemBinaryOutput:
    """
    Reads binary Q-Chem files from scratch or save directory.

    If Q-Chem is run with the @c -save option, the binary scratch files that are
    produced will be retained in a specified subdirectory. The data in them can
    be parsed as raw binary data.
    """

    # Enumerated file numbers in Q-Chem
    #
    # These are NOT documented officially!

    ## File number for (sparse) overlap derivative matrix
    FILE_SX      = 23
    ## File number for density matrices (both alpha and beta)
    FILE_DENSITY_MATRIX=54
    ## File number for overlap matrix
    FILE_OVERLAP_MATRIX=320

    def __init__(self, path = ''):
        """
        @param path Relative path to Q-Chem save directory
        """
        self.path = path

    def OverlapMatrix(self, n = None):
        """
        Reads the overlap matrix.
        
        @param n Dimension of (square) matrix. Default: None, will guess based
        on file contents.
        """
        filename = os.path.join(self.path, str(self.FILE_OVERLAP_MATRIX)+'.0')
        M = numpy.fromfile(filename, numpy.float64)
        if n is None:
            n = (M.size)**0.5
        M = M.reshape((n,n))
        return M
    
    def DensityMatrix(self, SeparateSpin = False, n = None):
        """
        Reads the one-electron density matrix.
   
        @param SeparateSpin returns separate alpha and beta density matrices,
        otherwise returns their sum. Default: False

        @param n Dimension of (square) matrix. Default: None, will guess based
        on file contents.
        """
        filename = os.path.join(self.path, str(self.FILE_DENSITY_MATRIX)+'.0')
        M = numpy.fromfile(filename, numpy.float64)
        if n is None:
            n = (M.size/2)**0.5
        DensityMat = M.reshape((2, n, n))
        if SeparateSpin:
            return DensityMat[0,:,:], DensityMat[1,:,:]
        else:
            return DensityMat[0,:,:] + DensityMat[1,:,:]

    def OverlapMatrixDerivative(self):
        """
        Reads the overlap derivative matrix.
        @warning Currently does nothing.
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

