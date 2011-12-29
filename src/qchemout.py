#!/usr/bin/env python
"""
Low-level routines to deal with Q-Chem input decks and
output files

Jiahao Chen, 2009-12
non-update 2011-03
"""

_qchemcmd = 'qchem'

####################
# Code starts here #
####################
from numpy import array, hstack
import logging

logging.basicConfig()
logger = logging.getLogger(__name__)

class QChemInformation:
    def __init__(self):
        self.Version = None

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

            Handlers = [handler_FatalError(),
                handler_Input(),
                handler_Geometry(),
                handler_GaussianBasis(),
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
                            self.thisData = None
                            state = None
                            #print response[0]
                            #print response[1]
                    if state is None:
                        for registered_handler in Handlers:
                            if registered_handler.trigger(line):
                                state = registered_handler
                                break
                    if state is None: print state, line,

                #Flush the last handler
                if state is not None:
                    self.Data.append(state.flush())

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
                t = line.split()
                self.data.append((int(t[0]), float(t[1]), float(t[2])))
                idx += 5
                message = line[idx:].strip()
                if message != '': #Have message
                    if message != 'Convergence criterion met':
                        print 'ERROR', message


class handler_NormalTermination(handler):
    def trigger(self, line):
        return '*** MISSION COMPLETED -- STARFLEET OUT ***' in line
    def handler(self, line):
        self.data = True
        return self.flush()


class handleratomiccharges(handler):
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


# Run as standalone
if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        print QChemOutput(sys.argv[1])
