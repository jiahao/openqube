#!/usr/bin/env python
"""
Q-Chem text output parser

Jiahao Chen <jiahao@mit.edu> 2011-12-28
"""

####################
# Code starts here #
####################
from numpy import array, hstack, zeros
from numpy.linalg import norm
import logging

from data import eV

logging.basicConfig()
logger = logging.getLogger(__name__)

class QChemOutput:
    """
    Container for parsed data from a Q-Chem text output file.
    """
    def __init__(self, filename = '', doAutoParse = True):
        """
        Initialize by saving filename and checking that it exists.
        
        Immediately begins parsing file, unless otherwise specified.

        @param filename Name of Q-Chem output file. Default: ''
        @param doAutoParse whether to immediately begin parsing. Default: True
        @throws IOError if file does not exist when parse.
        """
        self.filename = filename
        self.thisData = None
        self.Data = []
        if doAutoParse:
            self.Parse()

    def __repr__(self):
        buf = []
        for name, data in self.Data:
            buf.append('Perceived data: '+name)
            buf.append(str(data))
        return '\n'.join(buf)

    def Parse(self, Handlers = None):
        """Finite state machine.
        
        @param Handlers list of handlers. Default: None, in which the
        internal handlers will be used.
        All handlers beginning with handler_*() will be used.
        """

        if Handlers is None:

            Handlers = [handler_QChemVersion(),
                handler_FatalError(),
                handler_Input(),
                handler_Geometry(),
                handler_GaussianBasis(),
                handler_ModelChemistry(),
                handler_NumberOfElectrons(),
                handler_CoreHamiltonianMatrix(),
                handler_MultipoleMatrix(),
                handler_OverlapMatrix(),
                handler_OrthonormalizationMatrix(),
                handler_KineticEnergyMatrix(),
                handler_NuclearAttractionMatrix(),
                handler_SCFConvergence(),
                handler_FinalAlphaMOEigenvalues(),
                handler_FinalBetaMOEigenvalues(),
                handler_FinalAlphaMOCoefficients(),
                handler_FinalBetaMOCoefficients(),
                handler_FinalAlphaDensityMatrix(),
                handler_FinalBetaDensityMatrix(),
                handler_CISConvergence(),
                handler_TDAExcitationEnergy(),
                handler_TDDFTExcitationEnergy(),
                handler_FEDCoupling(),
                handler_SCFEnergyGradient(),
                handler_MullikenAtomicCharges(),
                handler_Molden(),
                handler_NormalTermination(),
                ]
            #TODO Dynamically register handlers?
        try:
            state = None
            with open(self.filename) as outfile:
                for line in outfile:
                    if state is not None:
                        response = state.handler(line)
                        if response is not None:
                            self.Data.append(response)
                            state.__init__() #Reset
                            self.thisData = None
                            state = None
                    if state is None:
                        for registered_handler in Handlers:
                            if registered_handler.trigger(line):
                                state = registered_handler
                                break
                    #if state is None: print state, line,

                #Flush the last handler
                if state is not None:
                    self.Data.append(state.flush())
                    state.__init__() #Reset

        except IOError:
            error_msg = "Filename "+self.filename+" does not appear to be valid"
            logger.error(error_msg)
            raise IOError, error_msg

    
class handler:
    def __init__(self):
        self.data = None
    def trigger(self, line):
        return line != ''
    def handler(self, line):
        return line
    def flush(self):
        return self.__class__.__name__.replace('handler_',''), self.data

class handler_QChemVersion(handler):
    def trigger(self, line):
        if 'Q-Chem, Version' in line:
            self.data = line.split()[2][:-1]
            return True
        else:
            return False
    def handler(self, line):
        return self.flush()

class handler_FatalError(handler):
    """Reads out part of Q-Chem output file with fatal error.
    Returns an exception that contains the error message and
    everything after it.
    If no error, returns empty string."""

    def trigger(self, line):
        if 'fatal error' in line:
            self.data = [line]
            return True
        else:
            return False

    def handler(self, line):
        self.data.append(line)
        
    def flush(self):
        error = ''.join(self.data)
        if error != '':
            raise ValueError, "Error in Q-Chem output file:\n" + error
        return 'FatalError', error

class handler_Input(handler):
    "Reads in Q-Chem input"
    def trigger(self, line):
        return 'User input:' in line
    def handler(self, line):
        if 'User input:' in line:
            return
        if '-'*62 in line:
            if self.data is None:
                self.data = []
            else:
                self.data = ''.join(self.data)
                return self.flush()
        else:
            self.data.append(line)

class handler_Geometry(handler):
    def trigger(self, line):
        return 'Standard Nuclear Orientation (Angstroms)' in line
    def handler(self, line):
        if self.data is None:
            if '-'*52 in line:
                self.data = []
        else: 
            if '-'*52 in line:
                #XXX Convert to recarray
                #This doesn't work?!
                #self.data = array(self.data, dtype=[('Element', str), ('Coordinates', float, 3)])
                return self.flush()
            else:
                t = line.split()
                self.data.append([t[1],]+map(float, t[2:5]))


class handler_NumberOfElectrons(handler):
    def trigger(self, line):
        t = line.split()
        if t[:2] == ['There', 'are'] and t[3:5] == ['alpha','and'] and \
                t[6:] == ['beta', 'electrons']:
            self.data = int(t[2]), int(t[5])
            return True
        else:
            return False

    def handler(self, line):
        return self.flush()


class handler_Molden(handler):
    """Reads a Molden input file from the Q-Chem output file"""
    def trigger(self, line):
        if '======= MOLDEN-FORMATTED INPUT FILE FOLLOWS =======' in line:
            self.data = []
            return True
        else:
            return False

    def handler(self, line):
        if '======= END OF MOLDEN-FORMATTED INPUT FILE =======' in line:
            return self.flush()

        self.data.append(line)
    
    def flush(self):
        return 'Molden', ''.join(self.data)


class handler_GaussianBasis(handler):
    """Counts the number of basis functions associated with atoms in
    AtomList as printed in the Q-Chem output file.

    If called with no AtomList, counts everything

    Populates self.BasisCount if not already populated.
    
    @todo Restore counting properties for atomlist
    @todo restore crash on PURECART"""

    def __init__(self):
        handler.__init__(self)
        self.BasisCount = None
        self.BasisName = ''
        self.dspher = True #Use spherical d Gaussians
        self.NumBasis = None
        self.NumBasis2= None
        self.thisAtom = None
        self.thisAtomBasisCount = None
        self.thisL = None

    def trigger(self, line):
        #if 'PURECART' in line.upper():
        #        error_msg = "Don't know how to handle PURECART"
        #        logger.error(error_msg)
        #        raise ValueError, error_msg

        if 'Requested basis set is' in line:
            self.BasisName = line.split()[-1]
            if '6-31' in self.BasisName and '6-311' not in self.BasisName:
                logger.debug("Detected Pople double zeta basis set. \
ing Cartesian d Gaussians")
                self.dspher = False
            return True

        else:
            return False

    def handler(self, line):
        if 'Total QAlloc Memory Limit' in line:
            return self.flush()

        elif 'There are' in line and 'basis functions' in line:
            self.NumBasis2 = int(line.split()[-3])

        elif 'Atom   I     L     Exponents     Normalized Contraction \
Coefficients' in line:
            self.BasisCount = []
        elif '-'*68 in line and self.thisAtom is not None:
            return self.flush()
        else:
            atomid = line[:5].strip()
            #The only other information we need is the angular momentum
            #quantum number
            L = line[12:19].strip()

            try: #Update if there is a new atom
                self.thisAtom = int(atomid)
                if self.thisAtomBasisCount is not None:
                    self.BasisCount.append(self.thisAtomBasisCount)
                self.thisAtomBasisCount = 0
            except ValueError:
                pass

            if L == '':
                pass
            elif L == '0-1':
                self.thisAtomBasisCount += 4
            else:
                try:
                    L = int(L)
                except ValueError:
                    return
               
                if L == 2: #d Gaussians
                    if self.dspher: #Using spherical d Gaussians?
                        self.thisAtomBasisCount += 5
                    else:
                        self.thisAtomBasisCount += 6
                elif L < 5:
                    #Q-Chem uses spherical Gaussians for d,f,g
                    self.thisAtomBasisCount += 2 * L + 1
                elif L >= 5:
                    #Q-Chem uses Cartesian Gaussians for h
                    self.thisAtomBasisCount += (L + 1) * (L + 2) / 2

    def flush(self):
        if self.thisAtom is None: #No explicit basis in output file
            return 'Basis', self.BasisName
        
        else:
            self.BasisCount.append(self.thisAtomBasisCount)
            assert self.NumBasis2 == sum(self.BasisCount)
            return 'Basis', (self.BasisName, self.BasisCount)

        
        ##Done parsing, now calculate
        #if AtomList == None:
        #    Total = sum(self.BasisCount)
        #else:
        #    Total = sum([self.BasisCount[x - 1] for x in AtomList])
        #return Total

class handler_ModelChemistry(handler):
    """
    Sample output:
        Exchange:     0.2500 Hartree-Fock + 0.7500 PBE
        Correlation:  1.0000 PBE
    """
    def trigger(self, line):
        if 'Exchange:' in line:
            self.data = {'Exchange':[], 'Correlation':[]}
            for x in line.split(':')[1].split('+'):
                x = x.strip()
                idx = x.find(' ')
                self.data['Exchange'].append((float(x[:idx]), x[idx+1:]))
            return True
    def handler(self, line):
        if 'Correlation:' in line:
            for x in line.split(':')[1].split('+'):
                x = x.strip()
                idx = x.find(' ')
                self.data['Correlation'].append((float(x[:idx]), x[idx+1:]))
        else:
            return self.flush()

class handler_SCFConvergence(handler):
    def trigger(self, line):
        return 'Cycle       Energy         DIIS Error' in line
    def handler(self, line):
        if '-'*39 in line:
            if self.data is None:
                self.data = []
            else:
                return self.flush()
        else:
            idx = line.find('00000')
            if idx > -1: #Has convergence data
                try:
                    t = line.split()
                    self.data.append((int(t[0]), float(t[1]), float(t[2])))
                    idx += 5
                    message = line[idx:].strip()
                    if message != '': #Have message
                        if message != 'Convergence criterion met':
                            print 'ERROR', message
                except ValueError:
                    pass

class handler_CISConvergence(handler):
    """
    CIS

    also called by TDDFT
    """
    def trigger(self, line):
        return 'Iter    Rts Conv    Rts Left    Ttl Dev     Max Dev' in line

    def handler(self, line):
        if '-'*51 in line:
            if self.data is None:
                self.data = []
            else:
                return self.flush()
        else:
            try:
                t = line.split()
                self.data.append((int(t[0]), int(t[1]), int(t[2]), float(t[3]),
                    float(t[4])))
                message = line[52:].strip()
                if message != '': #Have message
                    if message != 'Roots Converged':
                        print 'ERROR', message
            except (IndexError, ValueError):
                pass

class handler_NormalTermination(handler):
    def trigger(self, line):
        return '*** MISSION COMPLETED -- STARFLEET OUT ***' in line
    def handler(self, line):
        self.data = True
        return self.flush()

class ElectronicState:
    def __init__(self):
        self.Amplitudes = {} # Ideally we want this to be a sparse matrix!
                             # But we will make do with a dictionary[(OCC, VIRT)] = ampl
        self.Energy = None
        self.ExcitationEnergy = None
        self.Multiplicity = None
        self.Index = None
        self.Strength = None #Oscillator strength
        self.TransitionDipole = None #Transition dipole

    def __repr__(self):
        buf = [str(a)+' '+str(b) for a, b in self.__dict__.items()]
        return '\n'.join(buf)

    def isValid(self):
        "Run some sanity checks"
        #Test for correct unit of transition dipole
        #The oscillator strength is related to the energy and transition dipole (in a.u.) by
        #f = 2/3 * delta_E * trans_dipole**2
        #The claim is that we store everything in atomic units
        #so this can be used as a test to see if we got the right unit for the transition dipole
        if self.Strength > 0:
            test = 2.0/3 * self.ExcitationEnergy * norm(self.TransitionDipole)**2 / self.Strength
            assert abs(test - 1) < 0.15, 'Incorrect unit for transition dipole, discrepancy factor is '+str(test**0.5)
            
        return True

def dict_to_matrix(mydict, offset = 1):
    """Accepts a coordinate indexed dictionary
    and returns the corresponding dense matrix"""

    max_m, max_n = 1, 1
    for m, n in mydict.keys():
        max_m = max(m, max_m)
        max_n = max(n, max_n)
    M = zeros((max_m, max_n))
    for (m, n), x in mydict.items():
        M[m-offset, n-offset] = x
    return M
        

class handler_FEDCoupling(handler):
    """
          States        X12(D)      X12(A)        dX12    Coupling(eV)' 
        ------------------------------------------------------------------------------
          1    2     -0.000001    0.000002   -0.000003   -0.3642343E-06
          ...
    """
    def trigger(self, line):
        return 'States        X12(D)      X12(A)        dX12    Coupling(eV)' \
            in line

    def handler(self, line):
        if '-'*78 in line:
            if self.data is None:
                self.data = {}
            else:
                self.data = dict_to_matrix(self.data)
                return self.flush()
        else:
            try:
                t = line.split()
                state1, state2 = int(t[0]) - 1, int(t[1]) -1
                self.data[(state1, state2)] = float(t[-1])*eV
            except (IndexError, ValueError):
                pass

SpinMultiplicity = { 
    "Singlet":0,
    "Triplet":2,
}

class handlerexcitations(handler):
    """
    Parse entry of the form
        Excited state   2: excitation energy (eV) =    2.0813
           Total energy for state   2:  -3332.812081058783
           Trans. Mom.:  0.0371 X   0.0059 Y  -0.0349 Z
           Multiplicity: Singlet
            Strength   :  0.8578
           X: D(265) --> V(  1) amplitude = -0.2404
                 ...
    """
    def __init__(self, ExcitationName = None):
        handler.__init__(self)
        self.ExcitationName = ExcitationName
        self.State = ElectronicState()
    def trigger(self, line):
        return self.ExcitationName in line
    def handler(self, line):
        t = line.split()
        if 'Excited state' in line and 'excitation energy (eV)' in line:
            self.State = ElectronicState()
            self.State.Index, self.State.ExcitationEnergy = \
                int(t[2][:-1]), float(t[-1])*eV
        elif 'Total energy for state' in line:
            self.State.Energy = float(t[-1])
        elif 'Trans. Mom.:' in line:
            self.State.TransitionDipole = array(map(float, (t[2], t[4], t[6]))) 
        elif len(t)>0 and t[0] == 'Multiplicity':
            self.State.Multiplicity = SpinMultiplicity[t[-1]]
        elif len(t)>0 and t[0] == 'Strength':
            self.State.Strength = float(t[-1])
        elif '-->' in line:
            amplitude = float(t[-1])
            x = line.find('(')
            y = line.find(')')
            mo_occ = int(line[x+1:y])
            x = line.find('(',y+1)
            y = line.find(')',y+1)
            mo_virt = int(line[x+1:y])
            self.State.Amplitudes[mo_occ, mo_virt] = amplitude

        elif len(t) == 0:
            # Terminate
            assert self.State.isValid(), 'Error parsing electronic state'
            self.data.append(self.State)
        elif '-'*51 in line:
            if self.data is None:
                self.data = []
            else:
                return self.flush()

class handler_TDAExcitationEnergy(handlerexcitations):
    def __init__(self):
        handlerexcitations.__init__(self)
        self.ExcitationName = 'TDDFT/TDA Excitation Energies'

class handler_TDDFTExcitationEnergy(handlerexcitations):
    def __init__(self):
        handlerexcitations.__init__(self)
        self.ExcitationName = 'TDDFT Excitation Energies'

class handleratomiccharges(handler):
    def __init__(self, ChargeName = None):
        handler.__init__(self)
        self.ChargeName = ChargeName
    def trigger(self, line):
        return self.ChargeName in line
    def handler(self, line):
        if '-'*40 in line:
            if self.data is None:
                self.data = []
            else:
                self.data = array(self.data, dtype='a2,float64')
                return self.flush()
        else:
            try:
                t = line.split()
                self.data.append((t[1], float(t[2])))
            except (AttributeError, ValueError):
                pass

class handler_MullikenAtomicCharges(handleratomiccharges):
    def __init__(self):
        handleratomiccharges.__init__(self)
        self.ChargeName = 'Ground-State Mulliken Net Atomic Charges'

class handlermatrix(handler):
    """

    Sample output:
        Final Alpha density matrix.
                   1           2           3           4    
           1   1.0689856  -0.2711237   0.0489085  -0.0127631
           2  -0.2711237   1.0880717  -0.2101739   0.0187229
    ...
    """

    def __init__(self):
        handler.__init__(self)
        self.Buffer = None
        self.Matrix = None
        self.MatrixName = 'The Matrix'

    def trigger(self, line):
        return self.MatrixName in line

    def handler(self, line):
        """Parses a Q-Chem output file for a matrix."""
        w = line.split()
        #Are we done? Check if first and last tokens in line are numbers
        try:
            int(w[0])
            float(w[-1])
        except (IndexError, ValueError):
            #Done reading in stuff, this is other output
            #Transfer columns that have already been read in
            if self.Buffer is not None:
                if self.Matrix is None:
                    self.Matrix = array(self.Buffer)
                else:
                    self.Matrix = hstack((self.Matrix, array(self.Buffer)))
            return self.flush()

        #Check if it's a header or not
        #See if last entry contains a decimal point
        if '.' not in w[-1]:
            #Transfer columns that have already been read in
            if self.Buffer is not None:
                if self.Matrix is None:
                    self.Matrix = array(self.Buffer)
                else:
                    self.Matrix = hstack((self.Matrix, array(self.Buffer)))
            self.Buffer = []
        else: #This is not a header
            self.Buffer.append(map(float, w[1:]))

        

    def flush(self):
        return self.MatrixName, self.Matrix


class handler_CoreHamiltonianMatrix(handlermatrix):
    def __init__(self):
        handlermatrix.__init__(self)
        self.MatrixName = 'Core Hamiltonian Matrix'

class handler_MultipoleMatrix(handlermatrix):
    #XXX Multipole matrix order not correctly preserved since this is only a
    # partial match for the full name printed
    def __init__(self):
        handlermatrix.__init__(self)
        self.MatrixName = 'Multipole Matrix'

class handler_OverlapMatrix(handlermatrix):
    def __init__(self):
        handlermatrix.__init__(self)
        self.MatrixName = 'Overlap Matrix'

class handler_OrthonormalizationMatrix(handlermatrix):
    def __init__(self):
        handlermatrix.__init__(self)
        self.MatrixName = 'Orthonormalization Matrix'

class handler_KineticEnergyMatrix(handlermatrix):
    def __init__(self):
        handlermatrix.__init__(self)
        self.MatrixName = 'Kinetic Energy Matrix'

class handler_NuclearAttractionMatrix(handlermatrix):
    def __init__(self):
        handlermatrix.__init__(self)
        self.MatrixName = 'Nuclear Attraction Matrix'

class handler_FinalAlphaMOEigenvalues(handlermatrix):
    def __init__(self):
        handlermatrix.__init__(self)
        self.MatrixName = 'Final Alpha MO Eigenvalues'

class handler_FinalBetaMOEigenvalues(handlermatrix):
    def __init__(self):
        handlermatrix.__init__(self)
        self.MatrixName = 'Final BetaMO Eigenvalues'

class handler_FinalAlphaMOCoefficients(handlermatrix):
    def __init__(self):
        handlermatrix.__init__(self)
        self.MatrixName = 'Final Alpha MO Coefficients'

class handler_FinalBetaMOCoefficients(handlermatrix):
    def __init__(self):
        handlermatrix.__init__(self)
        self.MatrixName = 'Final Beta MO Coefficients'

class handler_FinalAlphaDensityMatrix(handlermatrix):
    """
    Reads a density matrix out of the output"
    
    @note Remember to put SCF_GUESS_PRINT 2 in Q-Chem $REM block
    @note reconstruct density matrix
    """
    def __init__(self):
        handlermatrix.__init__(self)
        self.MatrixName = 'Final Alpha density matrix'

class handler_FinalBetaDensityMatrix(handlermatrix):
    def __init__(self):
        handlermatrix.__init__(self)
        self.MatrixName = 'Final Beta density matrix'

class handler_SCFEnergyGradient(handlermatrix):
    def __init__(self):
        handlermatrix.__init__(self)
        self.MatrixName = 'Gradient of SCF Energy'


# Run as standalone
if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        print QChemOutput(sys.argv[1])
    else:
        print __doc__
        print 'To run on a particular Q-Chem output file, specify the filename'

