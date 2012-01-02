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

from data import Angstrom, eV

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
        self.Data = []
        try:
            open(self.filename)
            if doAutoParse:
                self.Parse()
        except IOError: #Must parse manually
            pass

    def __repr__(self):
        if self.filename == '':
            return 'No Q-Chem output loaded.'

        buf = ['Filename: '+self.filename]

        if len(self.Data) == 0:
            buf.append('\nNo Q-Chem output perceived.')

        for name, data in self.Data:
            buf.append('Perceived data: '+name)
            buf.append(str(data))

        return '\n'.join(buf)

    def Parse(self, Handlers = None):
        """Finite state machine.
        
        @param Handlers list of handlers. Default: None, in which the
        internal handlers will be used.
        All handler classes beginning with _handler_*() will be dynamically
        instantiated when parsed.
        """

        if Handlers is None:
            #Dynamically register handlers
            Handlers = [handler() for name, handler \
                    in globals().items() if name[:8]=='_handler']
        try:
            state = None
            with open(self.filename) as outfile:
                for line in outfile:
                    if state is not None:
                        response = state.handler(line)
                        if response is not None:
                            self.Data.append(response)
                            state = None
                    #If current handler finishes, recycles current line to
                    #see if it triggers another handler
                    if state is None:
                        for handler in Handlers:
                            if handler.trigger(line):
                                logger.debug('Triggered handler:', str(handler))
                                state = handler
                                break
                    #if state is None: print state, line,

                #Flush the last handler
                if state is not None:
                    self.Data.append(state.flush())

        except IOError:
            error_msg = "Filename "+self.filename+" does not appear to be valid"
            logger.error(error_msg)
            raise IOError, error_msg

class _superhandler:
    """
    Base class for Q-Chem output parser handler classes.
    """
    def __init__(self):
        ## Data to be saved
        self.data = None

    def trigger(self, line):
        """
        Criterion for triggering handler(). Oftentimes this will be aiming to
        detect the header of a table or something similar.
        
        @param line the current line of text to test
        @returns Boolean to see if anything in the current line is recognized
        as the start of some data block.
        """
        return line != ''

    def handler(self, line):
        """
        Parses a line of text.
        
        @param line the current line of text to parse
        @returns None each time it is called until it is ready to return data
        produced with flush()
        """
        if line != '':
            return self.flush()

    def flush(self):
        """
        Returns parsed data and resets internal data store.

        @returns a duple containing the name of the data being returned in the
        [0] element and the actual data in the [1] element.

        @note The default behavior is to return the name of the handler class
        with the beginning part '_handler_' removed, accompanied by the
        actual data. Unless otherwise specified, all derived handlers may be
        assumed to follow this convention.

        @note the flush operation is separated from the main handler in case
        the data needs to be flushed out manually.
        """
        data, self.data = self.data, None
        return self.__class__.__name__.replace('_handler_',''), data

class _handler_QChemVersion(_superhandler):
    """
    Parses the version of Q-Chem used to produce the output.
    
    @returns the version number as a string.

    Sample output parsed:
    @verbatim
 Q-Chem, Version 3.1, Q-Chem, Inc., Pittsburgh, PA (2007).
    @endverbatim
    """
    def trigger(self, line):
        if 'Q-Chem, Version' in line:
            self.data = line.split()[2][:-1]
            return True
        else:
            return False

    def handler(self, line):
        return self.flush()

class _handler_JobSeparator(_superhandler):
    """
    Parses the presence of a job separator in a batch job.
    
    @returns jobid if detected.

    Sample output parsed:
    @verbatim
*************************************************************
Job 2 of 2 
*************************************************************
    @endverbatim
    """
    def trigger(self, line):
        if 'Job ' in line and ' of ' in line:
            self.data = int(line.split()[1])
            return True
        else:
            return False

    def handler(self, line):
        return self.flush()

class _handler_FatalError(_superhandler):
    """
    Parses fatal errors in Q-Chem output.
    
    @returns The error message, or an empty string if there is no error.

    Sample output parsed:
    @verbatim
 Q-Chem fatal error occurred in module /usr/qchem/include/RefCount.h, line 48:

 RefCount::Constructor called with null pointer argument

 Sun Mar 20 22:01:19 2011
    @endverbatim
    """
    def trigger(self, line):
        if 'fatal error' in line:
            self.data = [line]
            return True
        else:
            return False

    def handler(self, line):
        self.data.append(line)
        
    def flush(self):
        self.data = ''.join(self.data)
        return _superhandler.flush(self)

class _handler_Input(_superhandler):
    """
    Parses Q-Chem input block echoed in output.
    
    @returns string containing the Q-Chem Input.

    @sa qchemin.py::QChemInput

    Sample output parsed:
    @verbatim
--------------------------------------------------------------
User input: 
--------------------------------------------------------------
$rem
exchange b3lyp
basis cc-pvtz
$end

$molecule
0 1
Na 0.000000 0.000000 -0.001000
Cl 0.000000 0.000000 2.169157
$end
    @endverbatim
    """
    def trigger(self, line):
        return 'User input:' in line
    def handler(self, line):
        if 'User input:' in line:
            return
        if '-'*62 in line:
            if self.data is None:
                self.data = []
            else:
                return self.flush()
        else:
            self.data.append(line)

    def flush(self):
        self.data = ''.join(self.data)
        return _superhandler.flush(self)

class _handler_Geometry(_superhandler):
    """
    Parses Q-Chem geometry in Cartesian coordinates
    
    @returns Cartesian geometry as a list of 4-tuples of the form
    ('Element', x, y, z)

    @todo this should really be returning a recarray in flush()
    Sample output parsed:
    @verbatim
 ----------------------------------------------------
       Standard Nuclear Orientation (Angstroms)
    I     Atom         X            Y            Z    
 ----------------------------------------------------
    1      Al      4.643205    -0.112268     0.418334
    2      O       3.344659     1.145007     0.851281
 ----------------------------------------------------
    @endverbatim
    """
    def trigger(self, line):
        if 'Standard Nuclear Orientation (Angstroms)' in line:
            self.Unit = Angstrom
            return True
        elif 'Standard Nuclear Orientation (Bohr)' in line:
            self.Unit = 1.0
            return True

    def handler(self, line):
        if self.data is None:
            if '-'*52 in line:
                self.data = []
        else: 
            if '-'*52 in line:
                #This conversion to recarray doesn't work?!
                #self.data = array(self.data, dtype=[('Element', str), 
                # ('Coordinates', float, 3)])
                return self.flush()
            else:
                t = line.split()
                self.data.append([t[1],]+ \
                        map(lambda x: float(x)*self.Unit, t[2:5]))


class _handler_NumberOfElectrons(_superhandler):
    """
    Parses the number of electrons of each spin.

    @returns Tuple of two integers being the number of alpha and beta
    electrons respectively.

    Sample output parsed:
    @verbatim
 There are       75 alpha and       75 beta electrons
    @endverbatim
    """
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


class _handler_Molden(_superhandler):
    """
    Parses a Molden input file from the Q-Chem output file.
    
    @returns the Molden input file as a string.

    Sample output parsed:
    @verbatim
======= MOLDEN-FORMATTED INPUT FILE FOLLOWS =======
[Molden Format]
[Atoms] (Angs)
 Al      1   13     10.93902031     -0.01715823     -0.75372654
  O      2    8      9.55711265      1.04162098     -1.40466878
  ...
======= END OF MOLDEN-FORMATTED INPUT FILE =======
    @endverbatim
    """
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
        self.data = ''.join(self.data)
        return _superhandler.flush(self)

class _handler_GaussianBasis(_superhandler):
    """
    Parses the basis set and, if present, number of basis functions associated
    with each atom.

    @returns a string containing the name of the basis set, and optionally, as
    the [1] component of a tuple, a list with the count of basis functions for
    each atom.

    @note This code may fail or give the wrong answer if PURECART is specified
    in the input. This is due to the builtin presumptions of Q-Chem's default
    assignments of spherical vs. Cartesian Gaussian basis functions.

    Sample output parsed:
    @verbatim
 Requested basis set is cc-pVTZ
 There are 120 shells and 340 basis functions
 Number of Cartesian basis functions =    390
 Maximum angular momentum            =      3
 Maximum degree of contraction       =     11
 Maximum number of momentum cases    =      1
 --------------------------------------------------------------------
 Atom   I     L     Exponents     Normalized Contraction Coefficients
 --------------------------------------------------------------------
   1    1     0    2.055000E+05   4.676786E-01
                   3.078000E+04   8.743786E-01
...
        2     0    2.055000E+05   1.207961E-03
                   3.078000E+04   1.922737E-03
                   7.006000E+03   3.279856E-03
...
 --------------------------------------------------------------------
    @endverbatim
    """

    def __init__(self):
        _superhandler.__init__(self)
        self.BasisCount = None
        self.BasisName = ''
        self.dspher = True #Use spherical d Gaussians
        self.NumBasis = None
        self.NumBasis2= None
        self.thisAtom = None
        self.thisAtomBasisCount = None
        self.thisL = None

    def trigger(self, line):
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

        elif 'Job number =' in line: #No more output
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
            self.data = self.BasisName
        else:
            self.BasisCount.append(self.thisAtomBasisCount)
            assert self.NumBasis2 == sum(self.BasisCount)
            self.data = (self.BasisName, self.BasisCount)

        return _superhandler.flush(self)


class _handler_ModelChemistry(_superhandler):
    """
    Parses the model chemistry.

    @returns a dictionary containing the keys 'Exchange' and 'Correlation'
    and items of lists of tuples of the form (coefficient, 'functional name')

    Sample output parsed:
    @verbatim
 Exchange:     0.2500 Hartree-Fock + 0.7500 PBE
 Correlation:  1.0000 PBE
    @endverbatim
    Sample output parsed:
    @verbatim
 Exchange:  PBE      Correlation:  PBE
    @endverbatim
    """
    def trigger(self, line):
        if 'Exchange:' in line:
            self.data = {'Exchange':[], 'Correlation':[]}
            self.handler(line)
            return True

    def handler(self, line):
        if 'Exchange:' in line:
            exchange_part = line[line.find('Exchange:')+10:line.find('Correlation:')]
            if '+' in exchange_part: #Combination specified
                for component in line.split('+'):
                    for coefficient, functional in component.split():
                        self.data['Exchange'].append(float(coefficient), functional.strip())
            else:
                self.data['Exchange'].append((1.0, exchange_part.strip()))

        if 'Correlation:' in line:
            line = line[line.find('Correlation:'):]
            if '+' in line: #Combination specified
                for component in line.split('+'):
                    for coefficient, functional in component.split():
                        self.data['Correlation'].append(float(coefficient), functional.strip())
            else:
                self.data['Correlation'].append((1.0, line.strip()))
        
        else:
            return self.flush()

class _handler_SCFConvergence(_superhandler):
    """
    Parses the SCF convergence history.

    @returns a list containing the cycle number, SCF energy and DIIS Error

    Sample output parsed:
    @verbatim
 ---------------------------------------
  Cycle       Energy         DIIS Error
 ---------------------------------------
    1      -4.1390916654      3.79E-09  00000 Convergence criterion met
 ---------------------------------------
    @endverbatim
    """
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

class _handler_CISConvergence(_superhandler):
    """
    Parses convergence history for configuration interaction singles (CIS)-type
    calculations.

    @returns a list of 5-tuples containing the iteration number, number of roots
    converged, number of roots remaining, total deviation and maximum deviation.

    @note such calculations include TDDFT
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
                    pass #Ignored for now
            except (IndexError, ValueError):
                pass

class _handler_NormalTermination(_superhandler):
    """
    Parses the existence of normal termination of a calculation

    @returns True

    Sample output parsed:
    @verbatim
        *************************************************************
        *                                                           *
        *  Thank you very much for using Q-Chem.  Have a nice day.  *
        *                                                           *
        *************************************************************


*** MISSION COMPLETED -- STARFLEET OUT ***
    @endverbatim    
    """
    def trigger(self, line):
        return '*** MISSION COMPLETED -- STARFLEET OUT ***' in line
    def handler(self, line):
        self.data = True
        return self.flush()


## Lookup table for spin quantum number of state given the English word.
SpinMultiplicity = { 
    "Singlet":0,
    "Doublet":1,
    "Triplet":2,
}

def dict_to_matrix(mydict, offset = 1):
    """
    Converts dictionary defining a coordinate-indexed sparse array and returns
    the corresponding dense array. Automatically creates the smallest possible
    array to hold the largest-indexed entry. Currently only works for
    two-dimensional arrays.

    @param mydict Dictionary to convert
    @param offset Whether the dictionary keys should be corrected for fencepost.
    Default: 0. Set to 1 if the coordinates are 1-indexed rather than 0-indexed.
    @returns numpy.ndarray
    """

    max_m, max_n = 1, 1
    for m, n in mydict.keys():
        max_m = max(m, max_m)
        max_n = max(n, max_n)
    M = zeros((max_m, max_n))
    for (m, n), x in mydict.items():
        M[m-offset, n-offset] = x
    return M
        

class _handler_FEDCoupling(_superhandler):
    """
    Parses couplings between electronic states as calculated using the Fragment
    Excitation Difference (FED) method.
   
    @sa http://pubs.acs.org/doi/abs/10.1021/jp076512i

    @returns a dictionary with each key being duple of integers @f$(m, n)@f$
    referring to the states being coupled and each corresponding value being a
    duple of @f$\Delta X_{mn} @f$ and the coupling matrix element (in a.u.)
    
    @f$\Delta X_{mn} @f$ is the difference in excitation number of the acceptor
    relative to the donor for the electronic transition @f$ m \rightarrow n @f$.
    The excitation number is a quantity with dimension of charge that represents
    the sum of the effective numbers of electrons and holes formed.

    Sample output parsed:
    @verbatim
            Fragment Excitations of Transition Densities and
              FED Couplings Between Singlet Excited States
 ------------------------------------------------------------------------------
   States        X12(D)      X12(A)        dX12    Coupling(eV)
 ------------------------------------------------------------------------------
   1    2      0.082846   -0.082846    0.165693    0.7941258E-04
 ------------------------------------------------------------------------------
    @endverbatim
    """
    def trigger(self, line):
        return 'States        X12(D)      X12(A)        dX12    Coupling(eV)' \
            in line

    def handler(self, line):
        if '-'*78 in line:
            if self.data is None:
                self.data = dict()
            else:
                #self.data = dict_to_matrix(self.data)
                return self.flush()
        else:
            try:
                t = line.split()
                state1, state2 = int(t[0]) - 1, int(t[1]) -1
                self.data[(state1, state2)] = float(t[-2]), float(t[-1])*eV
            except (IndexError, ValueError):
                pass


class ElectronicState:
    """
    Data structure containing information about electronic excited states.

    @note In spectroscopic terms, these are vertical excitations not adiabatic
    ones.
    """
    def __init__(self):
        ## Transition amplitudes stored as a dictionary with keys (occupied,
        #  virtual) and value being the transition amplitude.
        self.Amplitudes = dict()     
        ## Absolute energy @f$ E @f$
        self.Energy = None
        ## Excitation energy relative to the ground state, @f$ \Delta E @f$
        self.ExcitationEnergy = None
        ## Spin multiplicity @f$ S @f$ of the state. Stored as the spin quantum
        # number. @sa #SpinMultiplicity
        self.Multiplicity = None
        ## Index of state indicating relative order above the ground state
        self.Index = None
        ## Oscillator strength for the transtion from the ground state, @f$ f @f$
        self.OscillatorStrength = None
        ## Transition dipole from the ground state, @f$ \mu @f$
        self.TransitionDipole = None

    def __repr__(self):
        buf = [str(a)+' '+str(b) for a, b in self.__dict__.items()]
        return '\n'.join(buf)

    def isValid(self):
        """
        Run some sanity checks

        Currently, runs a simple test to see if the transition dipole has
        correct units.
        
        The oscillator strength is related to the energy and transition dipole
        (in a.u.) by the definition
        @f[
        f = \frac 2 3 \Delta E \vert \mu \vert^2
        @f]

        The claim is that we store everything in atomic units so this can be
        used as a test to see if we got the right unit for the transition
        dipole.
        """

        if self.OscillatorStrength > 0.0001:
            test = 2.0/3 * self.ExcitationEnergy * norm(self.TransitionDipole)**2 / self.OscillatorStrength
            assert abs(test - 1) < 0.15, 'Incorrect unit for transition dipole, discrepancy factor is '+str(test**0.5)
            
        return True


class _superhandler_ExcitedStates(_superhandler):
    """
    Superclass for parsers of excited state information.

    @returns list of ElectronicState instances

    Sample output parsed:
    @verbatim
 Excited state   1: excitation energy (eV) =    2.6020
    Total energy for state   1:  -3325.517023390698
    Multiplicity: Singlet
    Trans. Mom.:  0.0092 X   0.4489 Y  -0.0651 Z
    Strength   :  0.0131
    D(236) --> V(  2) amplitude =  0.8011
    D(238) --> V(  2) amplitude =  0.5895
    @endverbatim
    """
    def __init__(self, ExcitationName = None):
        """
        @param ExcitationName Electronic structure method used for the excited
        state calculation. Used to specify the text header that triggers the
        parser.
        """
        _superhandler.__init__(self)
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
            self.State.OscillatorStrength = float(t[-1])
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

class _handler_TDAExcitedStates(_superhandler_ExcitedStates):
    """
    Parses excited states calculated using the Tamm-Dancoff approximation (TDA)

    Sample output parsed:
    @verbatim
 ---------------------------------------------------
         TDDFT/TDA Excitation Energies
 ---------------------------------------------------

 Excited state   1: excitation energy (eV) =    2.6020
    Total energy for state   1:  -3325.517023390698
    Multiplicity: Singlet
    Trans. Mom.:  0.0092 X   0.4489 Y  -0.0651 Z
    Strength   :  0.0131
    D(236) --> V(  2) amplitude =  0.8011
    D(238) --> V(  2) amplitude =  0.5895

 ---------------------------------------------------
    @endverbatim
    """    
    def __init__(self):
        _superhandler_ExcitedStates.__init__(self)
        self.ExcitationName = 'TDDFT/TDA Excitation Energies'

class _handler_TDDFTExcitedStates(_superhandler_ExcitedStates):
    """
    Parses excited states calculated using time-dependent density functional
    theory (TDDFT)

    Sample output parsed:
    @verbatim
 ---------------------------------------------------
             TDDFT Excitation Energies
 ---------------------------------------------------

 Excited state   1: excitation energy (eV) =    2.6020
    Total energy for state   1:  -3325.517023390698
    Multiplicity: Singlet
    Trans. Mom.:  0.0092 X   0.4489 Y  -0.0651 Z
    Strength   :  0.0131
    D(236) --> V(  2) amplitude =  0.8011
    D(238) --> V(  2) amplitude =  0.5895

 ---------------------------------------------------
    @endverbatim
    """    
    def __init__(self):
        _superhandler_ExcitedStates.__init__(self)
        self.ExcitationName = 'TDDFT Excitation Energies'

class _superhandler_AtomicCharges(_superhandler):
    """
    Superclass for parsers of atomic charges.
    
    Sample output parsed:
    @verbatim
     Atom                 Charge (a.u.)
  ----------------------------------------
      1 Al                    1.279237
      2 O                    -1.279237
  ----------------------------------------
  Sum of atomic charges =     0.000000
    @endverbatim
    """
    def __init__(self, ChargeName = None):
        _superhandler.__init__(self)
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

class _handler_MullikenAtomicCharges(_superhandler_AtomicCharges):
    """
    Parses Mulliken atomic charges.
    
    Sample output parsed:
    @verbatim

          Ground-State Mulliken Net Atomic Charges

     Atom                 Charge (a.u.)
  ----------------------------------------
      1 Al                    1.279237
      2 O                    -1.279237
  ----------------------------------------
  Sum of atomic charges =     0.000000
    @endverbatim
    """
    def __init__(self):
        _superhandler_AtomicCharges.__init__(self)
        self.ChargeName = 'Ground-State Mulliken Net Atomic Charges'

class _superhandler_matrix(_superhandler):
    """
    Superclass for parsers of matrices.

    @returns matrix as a @code numpy.ndarray

    Sample output parsed:
    @verbatim
  Final Alpha density matrix.
          1           2           3           4    
  1   1.0689856  -0.2711237   0.0489085  -0.0127631
  2  -0.2711237   1.0880717  -0.2101739   0.0187229
...
    @endverbatim
    """

    def __init__(self):
        _superhandler.__init__(self)
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


class _handler_CoreHamiltonianMatrix(_superhandler_matrix):
    """
    Parses the core Hamiltonian matrix.
    """
    def __init__(self):
        _superhandler_matrix.__init__(self)
        self.MatrixName = 'Core Hamiltonian Matrix'

class _handler_MultipoleMatrix(_superhandler_matrix):
    """
    Parses a multipole matrix.
    
    @warning Multipole matrix order not correctly preserved since this is only a
    partial match for the full name printed
    """
    def __init__(self):
        _superhandler_matrix.__init__(self)
        self.MatrixName = 'Multipole Matrix'

class _handler_OverlapMatrix(_superhandler_matrix):
    """
    Parses an overlap matrix.
    """
    def __init__(self):
        _superhandler_matrix.__init__(self)
        self.MatrixName = 'Overlap Matrix'

class _handler_OrthonormalizationMatrix(_superhandler_matrix):
    """
    Parses an orthonormalization matrix.
    """
    def __init__(self):
        _superhandler_matrix.__init__(self)
        self.MatrixName = 'Orthonormalization Matrix'

class _handler_KineticEnergyMatrix(_superhandler_matrix):
    """
    Parses a kinetic energy matrix.
    """
    def __init__(self):
        _superhandler_matrix.__init__(self)
        self.MatrixName = 'Kinetic Energy Matrix'

class _handler_NuclearAttractionMatrix(_superhandler_matrix):
    """
    Parses a nuclear attraction matrix.
    """
    def __init__(self):
        _superhandler_matrix.__init__(self)
        self.MatrixName = 'Nuclear Attraction Matrix'

class _handler_FinalAlphaMOEigenvalues(_superhandler_matrix):
    """
    Parses converged eigenvalues corresponding to molecular orbitals of alpha
    spin.
    """
    def __init__(self):
        _superhandler_matrix.__init__(self)
        self.MatrixName = 'Final Alpha MO Eigenvalues'

class _handler_FinalBetaMOEigenvalues(_superhandler_matrix):
    """
    Parses converged eigenvalues corresponding to molecular orbitals of beta
    spin.

    @note This will not be printed for spin-restricted calculations.
    """
    def __init__(self):
        _superhandler_matrix.__init__(self)
        self.MatrixName = 'Final BetaMO Eigenvalues'

class _handler_FinalAlphaMOCoefficients(_superhandler_matrix):
    """
    Parses converged eigenvectors corresponding to molecular orbitals of alpha
    spin.
    """
    def __init__(self):
        _superhandler_matrix.__init__(self)
        self.MatrixName = 'Final Alpha MO Coefficients'

class _handler_FinalBetaMOCoefficients(_superhandler_matrix):
    """
    Parses converged eigenvectors corresponding to molecular orbitals of alpha
    spin.

    @note This will not be printed for spin-restricted calculations.
    """
    def __init__(self):
        _superhandler_matrix.__init__(self)
        self.MatrixName = 'Final Beta MO Coefficients'

class _handler_FinalAlphaDensityMatrix(_superhandler_matrix):
    """
    Parses a converged density matrix for electrons of alpha spin.
    
    @note The density matrix is printed only if @c SCF_GUESS_PRINT=2 is
    specified in the Q-Chem @c $REM block.
    """
    def __init__(self):
        _superhandler_matrix.__init__(self)
        self.MatrixName = 'Final Alpha density matrix'

class _handler_FinalBetaDensityMatrix(_superhandler_matrix):
    """
    Parses a converged density matrix for electrons of beta spin.
    
    @note This will not be printed for spin-restricted calculations.
    """
    def __init__(self):
        _superhandler_matrix.__init__(self)
        self.MatrixName = 'Final Beta density matrix'

class _handler_SCFEnergyGradient(_superhandler_matrix):
    """
    Parses an SCF Energy gradient.
    """
    def __init__(self):
        _superhandler_matrix.__init__(self)
        self.MatrixName = 'Gradient of SCF Energy'


# Run as standalone
if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        print QChemOutput(sys.argv[1])
    else:
        print __doc__
        print 'To run on a particular Q-Chem output file, specify the filename'

