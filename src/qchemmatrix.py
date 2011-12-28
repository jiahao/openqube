#!/usr/bin/env python
"""Q-Chem interface to calculate matrices

Jiahao Chen, 2009-12"""

from copy import deepcopy
import os, scipy

from data import Elements
from qchemout import QChemOutput
from qchemoutbin import QChemBinaryOutput

class FragmentData:
    "Container for nuclear and basis function information"
    def __init__(self, Fragment = None, BasisCount = None, FragZ = 0):
        """Initializes with Fragment, BasisCount and FragZ data"""
        if Fragment is None: self.Fragment = Fragment
        if BasisCount is None: self.BasisCount = BasisCount
        self.Z = FragZ



    def numFragBasisFns(self):
        """Returns the number of basis functions on the fragment."""
        N = 0
        for AtomId, ThisAtomNumBasis in enumerate(self.BasisCount):
            if (AtomId + 1) in self.Fragment:
                N += ThisAtomNumBasis
        return N



    def numTotalBasisFns(self):
        """Returns total number of basis functions in entire system"""
        return sum(self.BasisCount)



    def ProjectionMatrix(self):
        """Forms a projection matrix for the space spanned by the basis \
functions of the fragment.
       
Returns
-------
        
A square scipy.matrix containing the projection matrix.
"""
        n = sum(self.BasisCount)
        PI = scipy.matrix(scipy.zeros((n, n)))
        i = 0
        for AtomId, ThisAtomNumBasis in enumerate(self.BasisCount):
            for _ in range(ThisAtomNumBasis):
                if (AtomId + 1) in self.Fragment:
                    PI[i, i] = 1.
                i += 1
        return PI



    def CountFragZ(self, QChemData):
        """Counts total nuclear charge in fragment

Input
-----
QChemData - A QChemInputData object
"""
        FragZ = 0.
        for AtomId, atomline in enumerate(QChemData.read_input_section('molecule').strip().split('\n')[1:]):
            if (AtomId + 1) in self.Fragment:
                try:
                    Element = atomline.split()[0]
                    FragZ = int(Element)
                except ValueError:
                    try:
                        FragZ = Elements.index(Element.upper()) + 1
                    except ValueError:
                        print
                        raise ValueError, 'Unknown element symbol ' + Element
        self.Z = FragZ



    def NumAtomsTotal(self):
        """Counts total number of atoms in entire system"""
        return len(self.BasisCount)



    def FragmentBasisIterator(self):
        "Iterates over the basis function ids that are in the fragment"
        NowBasisID = 0
        for FragAtomID in self.Fragment:
            for thisAtomID, NumAtomBasis in enumerate(self.BasisCount):
                if FragAtomID == thisAtomID + 1:
                    for _ in range(NumAtomBasis):
                        yield NowBasisID
                        NowBasisID += 1
                else:
                    NowBasisID += NumAtomBasis

    def __iter__(self):
        yield self.FragmentBasisIterator()


def CalculateP(QChemData, FragmentList, QChemFileName = None, doOverwrite = True):
    """Gets Q-Chem to calculate the density matrix and overlap matrix
    and get them into a usable matrix within Python.
    
    Also populates information about Fragment.
    """

    #################################################
    # Make Q-Chem input deck, run it and parse output
    #################################################

    Input = deepcopy(QChemData)
    if QChemFileName == None:
        Input.filename = QChemData.filename[:-3] + '-stickytmp1.in'
    else:
        Input.filename = QChemFileName
    assert QChemData.filename != Input.filename, 'ERROR: Running this \
calculation would overwrite the original input file'

    OutputFilename = Input.filename[:-2] + 'out'

    doRunQChem = doOverwrite

    #First determine path to Q-Chem save directory
    if 'QCSCRATCH' in os.environ:
        savepath = os.environ['QCSCRATCH']
    else:
        print 'Warning, unknown scratch location.'
        savepath = ''

    #Our name for savedir
    savedir = 'stickytmp1'
    savepath = os.path.join(savepath, savedir)

    if not os.path.exists(savepath):
        doRunQChem = True

    if not os.path.exists(OutputFilename):
        doRunQChem = True

    #Either load saved data or run Q-Chem
    if doRunQChem:
        Output = Input.execute(savedir=savedir)
    else:
        Output = QChemOutput(OutputFilename)

    Output.CheckFatalError()

    n = sum(Output.BasisCount)

    QCBin = QChemBinaryOutput(savepath)
    S = QCBin.OverlapMatrix(n)
    P = QCBin.DensityMatrix(n)

    ##########################
    # Populate fragment info #
    ##########################

    Fragment = FragmentData()
    Fragment.Fragment = FragmentList
    Fragment.BasisCount = Output.BasisCount
    Fragment.CountFragZ(Input)

    ########
    # Done #
    ########
    return S, P, Fragment

