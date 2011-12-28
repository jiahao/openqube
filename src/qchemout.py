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
import logging, os, scipy

logging.basicConfig()
logger = logging.getLogger(__name__)

class QChemOutput:
    def __init__(self, filename = ''):
        #Check that file exists
        self.filename = filename

        if os.path.isfile(filename):
            #Count basis functions
            self.BasisCount = None
            self.CountBasisFunctions()
        else:
            error_msg = "Filename " + filename + " does not appear to be valid"
            logger.error(error_msg)
            raise IOError, error_msg


    def ReadMatrix(self, title, readrange = None):
        """Parses a Q-Chem output file for a matrix named _title_."""
        state = 'FindTitle'

        Matrix = []
        Buffer = []
        thisinstance = -1
        Matrices = []
        if readrange == None: readrange = [0]
        for line in open(self.filename):
            if state == 'FindTitle':
                if title in line:
                    thisinstance += 1
                    if thisinstance in readrange:
                        state = 'ReadCols'

            elif state == 'ReadCols':
                w = line.split()

                #Are we done? Check if first and last tokens in line are numbers
                try:
                    int(w[0])
                    float(w[-1])
                except (IndexError, ValueError):
                    #Done reading in stuff, this is other output

                    #Transfer columns that have already been read in
                    for MatrixCol in Buffer:
                        Matrix.append(MatrixCol)

                    #Done reading in Q-Chem output, now do sanity checking
                    Matrix = scipy.matrix(Matrix)
                    if Matrix.shape[1] == 0:
                        #Matrix not found
                        assert False, 'Error parsing instance '+\
                            str(thisinstance)+' of matrix '+title
                        Matrices.append(None)
                    else:
                        Matrices.append(Matrix)

                    if thisinstance == max(readrange): #Done
                        break
                    else: #Continue reading
                        state = 'FindTitle'

                #Check if it's a header or not
                #See if last entry contains a decimal point
                if '.' not in w[-1]: #This is a header

                    #Transfer columns that have already been read in
                    for MatrixCol in Buffer:
                        Matrix.append(MatrixCol)

                    #Read in new header
                    FirstColId = int(w[0])
                    LastColId = int(w[-1])

                    Buffer = []
                    for i in range(FirstColId, LastColId + 1):
                        Buffer.append([])

                else: #This is not a header

                    #Q-Chem's format is rownum __ __ ... __            

                    #Skip first entry
                    w.pop(0)

                    for i, entry in enumerate(w):
                        Buffer[i].append(float(entry))

        if len(Matrices) == 0:
            return
        elif len(Matrices) == 1:
            Matrices = Matrices[0]
        
        return Matrices


    def CheckFatalError(self):
        """Reads out part of Q-Chem output file with fatal error.
        Returns an exception that contains the error message and
        everything after it.
        If no error, returns empty string."""

        output = []

        accumulate = False
        for line in open(self.filename):
            if 'fatal error' in line:
                accumulate = True

            if accumulate:
                output.append(line)

        error = ''.join(output)
        if error != '':
            raise ValueError, "Error detected in Q-Chem output file " +\
                self.filename + '\n' + error



    def GetMolden(self):
        """Reads a Molden input file from the Q-Chem output file"""

        state = 0
        buf = []
        for line in open(self.filename):
            if '======= MOLDEN-FORMATTED INPUT FILE FOLLOWS =======' in line:
                state = 1
            elif '======= END OF MOLDEN-FORMATTED INPUT FILE =======' in line:
                break
            elif state == 1:
                buf.append(line)
        return ''.join(buf)



    def CountBasisFunctions(self, AtomList = None):
        """Counts the number of basis functions associated with atoms in
        AtomList as printed in the Q-Chem output file.

        If called with no AtomList, counts everything

        Populates self.BasisCount if not already populated."""

        if self.BasisCount == None:
            self.BasisCount = []

            dspher = True
            #True = use Spherical d Gaussians, else use Cartesian d Gaussians

            state = 'Find header'
            for line in open(self.filename):

                #Check for hard-coded use of Spherical vs
                #Cartesian Gaussians per basis set
                if 'Requested basis set' in line:
                    if '6-31' in line and '6-311' not in line:
                        logger.debug("Detected Pople double zeta basis set. \
Setting Cartesian d Gaussians")
                        dspher = False

                if state == 'Find header':
                    if "Atom   I     L     Exponents" in line:
                        state = 'Skip line'

                elif state == 'Skip line':
                    state = 'Count'
                    NumBasis = 0

                elif state == 'Count':
                    if '----------------' in line: #Read to the end
                        self.BasisCount.append(NumBasis)
                        self.BasisCount.pop(0)
                        break

                    atomid = line[:5].strip()
                    if atomid != '':
                        self.BasisCount.append(NumBasis)
                        NumBasis = 0

                    #The only other information we need is the angular momentum
                    #quantum number
                    L = line[12:19].strip()

                    if L == '0-1':
                        NumBasis += 4
                    elif L == '':
                        pass
                    else:
                        ll = int(L)
                        if ll == 2: #d Gaussians
                            if dspher: #Using spherical d Gaussians?
                                NumBasis += 5
                            else:
                                NumBasis += 6
                        elif ll < 5:
                            #Q-Chem uses spherical Gaussians for d,f,g
                            NumBasis += 2 * ll + 1
                        elif ll >= 5:
                            #Q-Chem uses Cartesian Gaussians for h
                            NumBasis += (ll + 1) * (ll + 2) / 2

                if 'PURECART' in line.upper():
                    error_msg = """\
Hey! You specified something for PURECART!"
I don't know what to do!"""
                    logger.error(error_msg)
                    raise ValueError, error_msg

        #Done parsing, now calculate
        if AtomList == None:
            Total = sum(self.BasisCount)
        else:
            Total = sum([self.BasisCount[x - 1] for x in AtomList])
        return Total



    def ReadDensityMatrix(self):
        "Reads a density matrix out of the output"

        ## For testing purposes, get initial guess matrices
        ## Remember to put SCF_GUESS_PRINT 2 in Q-Chem $REM block
        
        PA = self.ReadMatrix('Final Alpha density matrix')
        if PA is None:
            error_msg = "I don't know where the alpha density matrix is."
            logger.error(error_msg)
            raise ValueError, error_msg

        PB = self.ReadMatrix('Final Beta density matrix')
        if PB is None: #Can't find beta density matrix
            logger.debug("I don't know where the beta density matrix is; \
assuming it's the same as alpha.")
            PB = PA

        return PA + PB


###

def SanitizeRem(buf):
    newbuf = []
    keyword = ''
    for l in buf:
        if l[0] == '!':
            continue
        t = l.split()
        if len(t) >= 2:
            keyword = t[0]
            if keyword != 'jobtype':
                newbuf.append(l)

    return newbuf


def SanitizeMolecule(filename):
    """Sometimes Q-Chem fails to write correct optimized geometry block
    For now the sticky electron code only works with Cartesian formatted
    geometry
    Super cheat - use babel to generate correct geometry block"""

    os.system('babel -iqcout ' + filename + ' -oqcin TMP')
    buf = []
    state = None
    for l in open('TMP'):
        if '$end' in l:
            state = None
        elif '$molecule' in l:
            state = 'add'
        elif state == 'add':
            buf.append(l)
    return buf
