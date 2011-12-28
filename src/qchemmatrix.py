#!/usr/bin/env python
"""Q-Chem interface to calculate matrices

Jiahao Chen, 2009-12"""

from copy import deepcopy
import os
import scipy
import scipy.linalg as linalg
import itertools
import Stencil
from ChemicalData import Elements
from QChemIO import QChemOutput
from readbin import ReadOverlap, ReadDensity

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

    S = ReadOverlap(savepath, n)
    P = ReadDensity(savepath, n)

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



def CalculatePx(QChemData, Fragment, QChemFileName = None, \
    DistanceDisplaced = (0., 0., 0.0001), Verbose = True, DiffStencil = None):

    if DiffStencil is None:
        DiffStencil = Stencil.FirstOrderCentralDifferenceStencil()

    assert isinstance(DiffStencil, Stencil.Stencil)

    ##########################################
    #Step 1. Evaluate density matrices on grid
    ##########################################
    DensityMatricesOnGrid = []
    for n, _ in DiffStencil.coeffs:
        Input = deepcopy(QChemData)
        if QChemFileName is None:
            Input.filename = QChemData.filename[:-3] + '-stickytmp3-%d.in' % n
        else:
            Input.filename = QChemFileName
        assert QChemData.filename != Input.filename, 'ERROR: Running this calculation would overwrite the original input file'

        Input.write_input_section('molecule',
            PerturbFragment(Input.read_input_section('molecule'),
            Fragment, n * scipy.array(DistanceDisplaced)))

        Run = True
        OutputFilename = Input.filename[:-2] + 'out'
        if os.path.exists(OutputFilename):
            Output = QChemOutput(OutputFilename)
            try:
                Output.CheckFatalError()
                print 'Reusing existing output', OutputFilename
                DensityMatricesOnGrid.append(Output.ReadDensityMatrix())
                Run = False
            except ValueError:
                pass

        #First determine path to Q-Chem save directory
        if 'QCSCRATCH' in os.environ:
            savepath = os.environ['QCSCRATCH']
        else:
            print 'Warning, unknown scratch location.'
            savepath = ''

        #Our name for savedir
        savedir = 'stickytmp3'
        savepath = os.path.join(savepath, savedir)

        if Run:
            frags = QChemData.filename[:-3].split('-')
            Output = Input.execute(savedir=savedir)
            Output.CheckFatalError()

            DensityMatricesOnGrid.append(ReadDensity(savedir))

    ##################################
    #Step 2. Calculate density response
    ##################################
    Px = DiffStencil.ApplyToFunctionOnGrid(DensityMatricesOnGrid, \
            h = linalg.norm(DistanceDisplaced))

    ########
    # Done #
    ########
    return Px



def CalculateS(QChemData, DistanceDisplaced = (0, 0, 0.0001),
        QChemFileName = None, Verbose = True, DiffStencil = None):
    """Gets Q-Chem to calculate the overlap matrix for an artificially
    superposed dimerr and gets it into a usable matrix within Python.
    
    This function wraps various stencil approximants.
    """


    if DiffStencil is None:
        DiffStencil = Stencil.FirstOrderFivePointCentralDifferenceStencil()

    assert isinstance(DiffStencil, Stencil.Stencil)
    DiffStencil.AddZero()

    N = DiffStencil.GetNumPoints()

    ##################################
    #Step 1. Make Q-Chem input deck
    ##################################

    Input = deepcopy(QChemData)
    if QChemFileName is None:
        Input.filename = QChemData.filename[:-3] + '-stickytmp2.in'
    else:
        Input.filename = QChemFileName
    assert QChemData.filename != Input.filename, 'ERROR: Running this calculation would overwrite the original input file'

    #Make fake oligomer
    s = Input.read_input_section('molecule').strip()

    #Cook up an acceptable charge/spin multiplicity
    charge = N * int(s.split('\n')[0].split()[0])
    spin = int(s.split('\n')[0].split()[1])

    if N % 2 == 0: spin = 1

    Geometries = []
    All = range(1, len(s.split('\n'))) #A list of all the atoms
    for n, _ in DiffStencil.coeffs:
        Geometries.append(PerturbFragment(s, Fragment = All, \
            delta = n * scipy.asanyarray(DistanceDisplaced)).split('\n')[1:-1])

    Geometry = '%d %d\n' % (charge, spin) + '\n'.join(itertools.chain(*Geometries))
    Input.write_input_section('molecule', Geometry + '\n', Overwrite = True)
    Input.append_input_section('rem', "\nMAX_SCF_CYCLES 0 Don't solve SCF\n")

    ##################################
    #Step 2. Run Q-Chem
    ##################################

    #First determine path to Q-Chem save directory
    if 'QCLOCALSCR' in os.environ: #We expect this job to fail!
        del os.environ['QCLOCALSCR']
    if 'QCSCRATCH' in os.environ:
        savepath = os.environ['QCSCRATCH']
    else:
        print 'Warning, unknown scratch location.'
        savepath = ''

    #Our name for savedir
    savedir = 'stickytmp2'
    savepath = os.path.join(savepath, savedir)

    Run = True
    OutputFilename = Input.filename[:-2] + 'out'
    """
    if os.path.exists(OutputFilename):
        Output = QChemOutput(OutputFilename)
        try:
            #Output.CheckFatalError()
            print 'Reusing existing output', OutputFilename
            Run = False
            S = ReadOverlap(savepath)
        except ValueError:
            pass
    """

    if Run:
        frags = QChemData.filename[:-3].split('-')
        Output = Input.execute(savedir=savedir)

        ##################################
        #Step 3. Read overlap matrix
        ##################################

        S = ReadOverlap(savepath)

    ##################################
    #Step 4. Calculate the two types of overlap derivative
    ##################################

    #Need number of basis functions on monomer
    nc = S.shape[0] / N

    #Calculate derivative density matrices by finite difference
    Sx = scipy.matrix(scipy.zeros((nc, nc)))

    norm = linalg.norm(DistanceDisplaced)

    #Find zero
    iz = DiffStencil.GetIndex(0)

    for i, (_, c) in enumerate(DiffStencil.coeffs):
        Sx += c * S[(iz * nc):((iz + 1) * nc), (i * nc):((i + 1) * nc)]
    Sx /= norm

    #Calculate derivative density matrices by finite difference
    Sx0 = scipy.matrix(scipy.zeros((nc, nc)))
    for i, (_, c) in enumerate(DiffStencil.coeffs):
        Sx0 += c * S[(i * nc):((i + 1) * nc), (iz * nc):((iz + 1) * nc)]
    Sx0 /= norm

    #Sx is supposed to be antisymmetric
    #Purify by removing symmetric part
    if Verbose:
        print "Removed symmetric part of Sx, ratio =", \
        linalg.norm(Sx + Sx0) / linalg.norm(Sx - Sx0)

    Sx = 0.5 * (Sx - Sx0)

    #Form two-sided overlap derivative
    Sxx = scipy.matrix(scipy.zeros((nc, nc)))

    for i, (_, c) in enumerate(DiffStencil.coeffs):
        for j, (_, d) in enumerate(DiffStencil.coeffs):
                Sxx += c * d * S[(i * nc):((i + 1) * nc), (j * nc):((j + 1) * nc)]

    Sxx /= norm ** 2

    #Sxx is supposed to be symmetric
    #Purify by removing antisymmetric part
    if Verbose:
        print "Removed antisymmetric part of Sxx, ratio =", \
        linalg.norm(Sxx - Sxx.T) / linalg.norm(Sxx + Sxx.T)

    Sxx = 0.5 * (Sxx + Sxx.T)

    return Sx, Sxx



#######################################
# Functions for dealing with geometry
#

def PerturbFragment(geometry, Fragment = None, delta = None):
    if Fragment == None or delta == None:
        return geometry

    #Make displaced geometry
    buf = []
    atomid = 1
    for line in geometry.strip().split('\n'):
        #Does the line have three numbers after something else?
        try:
            w = line.split()
            x = float(w[1])
            y = float(w[2])
            z = float(w[3])
        except (IndexError, ValueError):
            buf.append(line)
        else:
            #Displace only the atoms in the fragment
            if atomid in Fragment:
                buf.append(w[0] + ' %f %f %f ' % (x + delta[0], y + delta[1], z + delta[2]))
            else:
                buf.append(line)

            atomid += 1

    return '\n'.join(buf) + '\n'
