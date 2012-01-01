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
import logging, os
from copy import deepcopy

logging.basicConfig()
logger = logging.getLogger(__name__)

class QChemInput:
    """
    """
    def __init__(self, InputDeckFilename = None):
        logger.debug('Initializing new QChemInput() with filename: %s',
                     InputDeckFilename)
        
        self.filename = InputDeckFilename
        self.jobs = []
        if InputDeckFilename is not None and \
                os.path.exists(InputDeckFilename):
            self.ImportQChemOutput(InputDeckFilename)


    def isValid(self):
        """
        Checks that data contained can be used to make a valid input deck

        @todo
        If the first job has READ in the geometry, make sure scratch file
        exists and -save is used.

        @sa Q-Chem 4.0 manual, section 3.2.1
        """
        #Delete possibly empty first job
        if repr(self.jobs[0]) == '':
            self.jobs.pop(0)
                
        return all([job.isValid() for job in self.jobs])



    def GetCurrentJob(self):
        if len(self.jobs) == 0: #Create empty job
            self.jobs.append(QChemInputJob())

        return self.jobs[-1]



    def write(self, InputDeckFilename = None):
        """
        Writes input file.

        @note The input deck will be necessary but possibly not sufficient to
        reconstruct the requested calculation exactly. Section 3.1 of the Q-Chem
        4.0 manual states that additional inputs absent from the input file may
        be taken from @c $HOME/.qchemrc, or @c $QC/config/preferences, or
        internal program defaults (which may be version-specific).

        @param InputDeckFilename Name of Q-Chem input file to write.
        Default: None; will return the entire output file as a string,
        i.e. same as calling __repr__()
        """

        if InputDeckFilename != None:
            self.filename = InputDeckFilename

        assert self.isValid(), "Insufficient data for valid Q-Chem input deck."

        buf = []
        for job in self.jobs:
            for input_section_keyword, input_section in job.sections:
                buf.append('$' + input_section_keyword + '\n')
                buf.append(input_section)
                buf.append('\n$end\n\n')
            buf.append('@@@\n')

        #Delete last end-of-job marker
        buf.pop(-1)

        outputtext = ''.join(buf)
        if self.filename == None:
            return outputtext
        else:
            with open(self.filename, 'w') as f:
                f.write(outputtext)


    def ImportQChemOutput(self, filename):
        job = QChemInputJob()
        job.filename = filename
        input_section_buffer = []
        input_section_keyword = None
        state = 'none'
        jobid = 0

        filenameroot = '.'.join(filename.split('.')[:-1])
        filenameext = filename.split('.')[-1]

        for line in open(filename):
            if '$end' in line:
                state = 'none'

                if input_section_keyword == 'molecule':
                    oldline1 = input_section_buffer[0]
                    #buf = SanitizeMolecule(filename)
                    #Q-Chem geometry optimization output is broken in 3.2:
                    #1. Final geometry is always output in some broken Z-matrix
                    #2. Number of electrons output is incorrect with ECP -__-
                    #manual override
                    input_section_buffer[0] = oldline1

                elif input_section_keyword == 'rem':
                    input_section_buffer = SanitizeRem(input_section_buffer)

                job.write_input_section(input_section_keyword, \
                     ''.join(input_section_buffer))

                input_section_buffer = []

            elif '$' in line:
                input_section_keyword = line.split()[0][1:]
                state = 'readblock'

            elif state == 'readblock':
                input_section_buffer.append(line)

            elif '@@@' in line:
                job.filename = '.'.join([filenameroot + '-job' + str(jobid),
                                         filenameext])
                jobid += 1
                self.jobs.append(job)
                job = QChemInputJob()

        job.filename = '.'.join([filenameroot + '-job' + str(jobid),
                                 filenameext])

        self.jobs.append(job)



    def execute(self, outfilename = None, qchemcmd = _qchemcmd,
                parse = False, savedir = None, Overwrite = False):
        """Runs Q-Chem input file. Returns a QChemOutput class."""

        #Commits current data to disk
        self.write()
        if outfilename is None:
            outfilename = self.filename[:-2] + 'out'

        if savedir is not None:
            #if 'QCLOCALSCR' in os.environ and 'QCSCRATCH' in os.environ \
            #        and os.environ['QCLOCALSCR'] == os.environ['QCSCRATCH']:
            #    del os.environ['QCLOCALSCR']

            dosave = '-save'
        else:
            dosave = ''
            savedir = ''


        OutputExists = os.path.exists(outfilename)
        if Overwrite or not OutputExists:
            from subprocess import Popen, STDOUT, PIPE
            try:
                cmdstring = ' '.join([qchemcmd, dosave, self.filename,
                    outfilename, savedir])
            except TypeError, e:
                logger.error("""Could not generate valid shell command line:
qchemcmd    = %s
filename    = %s
outfilename = %s
""", qchemcmd, self.filename, outfilename)
                raise e

            p = Popen(cmdstring, stdout = PIPE, stderr = STDOUT, shell = True)
            p.wait()
            print 'Console output is: %s' % p.stdout.read()
            logger.info('Console output is: %s', p.stdout.read())
        else:
            if OutputExists:
                logger.warning('Output file exists: %s and overwrite was not \
specified. Refusing to run Q-Chem.', outfilename)
        if parse:
            from qchemout import QChemOutput
            return QChemOutput(outfilename)



    def __repr__(self):
        return '\n@@@\n'.join([x.__repr__() for x in self.jobs])


    def __str__(self):
        return self.__repr__()

## Documented user input section keywords
# @sa Q-Chem 4.0 manual, Table 3.1 and Appendix C
KnownInputSectionKeywords = ["molecule", "rem", "end", "basis", "cdft", "comment",
"ecp", "empirical_dispersion", "external_charges", "force_field_params",
"intracule", "isotopes", "multipole_field", "nbo", "occupied", "opt", "pcm",
"pcm_solvent", "svp", "svpirf", "plots", "van_der_waals", "xc_functional"]

#rem keywords which can be read"
Reads = ["scf_guess", "geometry"]

class QChemRemArray:
    """
    Container class for the @c $rem array
    @sa Q-Chem 4.0 manual, Section 3.5
    
    Specifies type of calculation. with individual @c $rem variables.

    The format of this block consists of lines of the form
    @verbatim
    REM_VARIABLE VALUE [comment]
    @endverbatim

    @note Additional @c $rem variables may be read in from .qchemrc and
preferences.

    @note  A line beginning with @c ! is a comment
    """

class QChemInputJob:
    """

    In the geometry block:
    @todo
    Q-Chem ignores commas and equal signs, and requires all distances, positions and
    angles to be entered as Angstroms and degrees. unless @c INPUT_BOHR=TRUE in
    @c $rem , in which case all lengths are assumed to be in bohr.

    @sa Q-Chem 4.0 manual, Section 3.1.
    """

    #Here is a subclass for the Rem block
    """
    class QChemInputRemBlock(dict):
        def __init__(self, *args, **kwargs):

            if 'inputstring' in kwargs:
                inputstring = kwargs['inputstring']
                del kwargs['inputstring']
            else:
                inputstring = None

            dict.__init__(self, args, kwargs)
            
            if inputstring is None:
                return

            for line in inputstring.split('\n'):
                if line[0] == '!':
                    #Throw away comment
                    continue
                
                t = line.split()
                if len(t) >= 2:
                    keyword = t[0]
                    value = t[1]
                    dict.__setitem__(self, keyword, value)
                    #Throw away comments
    """
    def __init__(self, filename = None):
        self.sections = []
        self.filename = filename
        #self.rem = self.QChemInputRemBlock() #TODO finish migration to RemBlock

    def isValid(self):
        """
        Checks that data contained can be used to make a valid input deck
        @todo
        Section 3.2 of the Q-Chem 4.0 manual states that the net charge must be
        -50 <= Z <= 50 (0 for neutral) and the spin multiplicity must be
        1<=S<=10 (1 for a singlet).

        @todo Accept @c READ @c filename in the geometry block
        @sa Q-Chem 4.0 manual, Section 3.2.2

        @todo Check rules for Cartesian and Z-Matrix coordinates
        @sa Q-Chem 4.0 manual, Section 3.3-3.4


        """
        return self.has_input_section('rem') and \
               self.has_input_section('molecule')



    def has_input_section(self, input_section_keyword):
        return (input_section_keyword in [x for x, _ in self.sections])



    def read_input_section(self, input_section_keyword):
        for keyword, section in self.sections:
            if input_section_keyword == keyword:
                return section

        return None



    def delete_input_section(self, input_section_keyword):
        newsections = []
        for keyword, section in self.sections:
            if input_section_keyword != keyword:
                newsections.append([keyword, section])
        self.sections = newsections



    def write_input_section(self, input_section_keyword, input_section,
                            Overwrite = True):
        if Overwrite and self.has_input_section(input_section_keyword):
            self.delete_input_section(input_section_keyword)
        elif not Overwrite and self.has_input_section(input_section_keyword):
            return #Do nothing
        self.sections.append([input_section_keyword, input_section])

        #Special handling for REM block
        if input_section_keyword == 'rem':
            pass

    def append_input_section(self, input_section_keyword, input_section):

        #find requested input section keyword
        input_section_keywords = [x for x, _ in self.sections]
        n = -1
        for i, this_keyword in enumerate(input_section_keywords):
            if this_keyword == input_section_keyword:
                n = i
                break

        if n == -1: #If we got here, could not find input section
            logger.debug('Creating new input section: %s',
                         input_section_keyword)
            self.write_input_section(input_section_keyword, input_section)
        else:
            self.sections[n][1] += input_section

    def get_input_section(self, input_section_keyword):
        for a, b in self.sections:
            if a.lower().strip() == input_section_keyword.lower().strip():
                return b


    #Specific to REM section
    
    def rem_normalize(self):
        "Cleans up rem block"
        rem_section = self.get_input_section('rem')
        rem_dict = {}
        for l in rem_section.split('\n'):
            t = l.split()
            try:
                keyword, value, comments = t[0], t[1], ' '.join(t[2:]) 
                rem_dict[keyword] = value, comments
            except IndexError:
                pass

        maxlen = max([len(w) for w in rem_dict.keys()])
        maxlen2 = max([len(w) for w, _ in rem_dict.items()])

        new_rem_section = []

        for keyword, (value, comments) in rem_dict.items():
            padlen = maxlen - len(keyword) + 1
            padlen2 = maxlen2 - len(value) + 1
            new_rem_section.append(keyword + ' '*padlen + value + ' '*padlen2\
                                       + comments)

        self.write_input_section('rem', '\n'.join(new_rem_section), \
                              Overwrite = True)


    def rem_get(self, keyword):
        rem_section = self.get_input_section('rem')
        for l in rem_section.split('\n'):
            t = l.split()
            if len(t) >= 2:
                if keyword.upper() == t[0].upper():
                    return t[1]
        
        return None

    def rem_set(self, keyword, value):
        rem_section = self.get_input_section('rem')
        new_rem_section = []
        found = False
        for l in rem_section.split('\n'):
            t = l.split()
            if len(t) >= 2:
                thiskeyword, thisvalue = t[0], t[1]
                if len(t) > 2:
                    thiscomment = ' '.join(t[2:])
                else:
                    thiscomment = ''

                if keyword.upper() == thiskeyword.upper():
                    found = True
                    new_rem_section.append(' '.join((thiskeyword, value,
                                                    thiscomment)))
                else:
                    new_rem_section.append(' '.join((thiskeyword, thisvalue,
                                                    thiscomment)))
        
        if not found: #add it in
            new_rem_section.append(' '.join((keyword, value)))

        self.write_input_section('rem', '\n'.join(new_rem_section))


    def rem_delete(self, keyword):
        "Comments out keyword from rem block"
        rem_section = self.get_input_section('rem')
        new_rem_section = rem_section.replace(keyword, '!' + keyword)
        self.write_input_section('rem', new_rem_section, Overwrite = True)
        


    def execute(self, outfilename = None, qchemcmd = _qchemcmd, parse = True,
                savedir = '', Overwrite = False):
        """Runs Q-Chem input file. Returns a QChemOutput class."""

        if outfilename == None:
            outfilename = self.filename[:-2] + 'out'

        #Commits current data to disk
        tmpjob = QChemInput(self.filename)
        tmpjob.jobs.append(self)
        tmpjob.write()
        return tmpjob.execute(outfilename = outfilename, qchemcmd = qchemcmd,
                              Overwrite = Overwrite, savedir = savedir,
                              parse = parse)



    def write(self, filename = None, doCheck = True):
        """
        Writes input file
        
        filename name of filename to write.
        If None, returns as string. Same as __repr__()
        """

        if doCheck:
            assert self.isValid(), "Insufficient data for valid Q-Chem \
                    input deck."

        outputtext = self.__repr__()
        if filename == None:
            return outputtext
        else:
            with open(filename, 'w') as f:
                f.write(outputtext)

            self.filename = filename



    def __repr__(self):
        buf = []
        for input_section_keyword, input_section in self.sections:
            buf.append('$' + input_section_keyword + '\n')
            buf.append(input_section)
            buf.append('\n$end\n\n')

        return ''.join(buf)

def QChemInputForElectrostaticEmbedding(ResID, CHARMM_CARD_file,
    CHARMM_RTFile, InputDeck = None, QChemInputFileName = None):
    '''
    Create a Q-Chem input deck to calculate the electrostatic field at each MM
    atom due to the quantum-mechanical region.

    @param ResID CHARMM ResID which defines the embedded QM region 

    @param CHARMM_CARD_file Name of CHARMM CARD file containing molecular
    coordinates

    @param CHARMM_RTFile Name of CHARMM Residue Topology File containing
    force field parameters (notably, atomic charges) (string or iterable)

    @param InputDeck Q-Chem input deck. (Optional.) \
    If value is a string, is treated as a filename for QChemInput object. \
    If it is a QChemInput object, then it will be modified. \
    If not specified, a new QChemInput object will be created.

    @param NewFileName A new filename for InputDeck. (Optional.)

    @returns A Q-Chem input deck that is ready to be written and executed.

    NOTE: The file will NOT be written to disk. issue a .write() manually as
    necessary!               
    '''

    #Load fixed charge specification in RTF file
    ChargeParameter = {}
    try:
        for l in open(CHARMM_RTFile):
            if l[:4].upper() == 'ATOM':
                t = l.split()
                AtomType, Charge = t[1].upper(), float(t[3])
                ChargeParameter[AtomType] = Charge

    except TypeError:
        #Not a string, assume it's an iterable of strings
        for filename in CHARMM_RTFile:
            for l in open(filename):
                if l[:4].upper() == 'ATOM':
                    t = l.split()
                    AtomType, Charge = t[1].upper(), float(t[3])
                    ChargeParameter[AtomType] = Charge

    #Generate $external_charges and $molecule blocks for Q-Chem input
    QMbuf = ['0 1'] #XXX Hard-coded charge and spin!
    MMbuf = []
    for l in open(CHARMM_CARD_file):
        t = l.split()
        if len(t) > 8: #Need at least 9 columns
            AtomType, x, y, z, thisResID = t[3].upper(), float(t[4]), \
                float(t[5]), float(t[6]), t[8]
            if thisResID == ResID: #QM region
                QMbuf.append(AtomType[0]+('%15.8f'*3) % (x, y, z))
            else:
                assert AtomType in ChargeParameter, 'ERROR: Could not find \
charge parameter for '+AtomType+' in '+CHARMM_RTFile
                charge = ChargeParameter[AtomType]
                MMbuf.append(('%15.8f\t'*4) % (x, y, z, charge)) 

    #Make Q-Chem input deck
    if InputDeck == None:
        InputDeck = QChemInput(QChemInputFileName)
    elif type(InputDeck) == type('1'): #A string, treat as filename
        InputDeck = QChemInput(InputDeck)
    elif isinstance(InputDeck, QChemInput):
        pass
    else:
        assert False, "I don't know what to do with InputDeck ="+str(InputDeck)
    InputDeck.GetCurrentJob().append_input_section('rem', 'igdesp\t%d\n' %
                                                   len(MMbuf))
    InputDeck.GetCurrentJob().write_input_section('molecule', '\n'.join(QMbuf),
                                                  Overwrite = True)
    InputDeck.GetCurrentJob().write_input_section('external_charges',
                                    '\n'.join(MMbuf), Overwrite = True)
    if QChemInputFileName != None: InputDeck.filename = QChemInputFileName
    return InputDeck #Warning, not written to disk yet!!



def QChemInputForTransitionDipole(InputDeck1, InputDeck2,
                                  NewFileName = 'tmp.in'):
    '''
    Create a Q-Chem input deck to calculate the transition dipole
    moment between two electronic states.

    This makes use of CDFT-CI's ability to calculate transition dipole moments.

    @param InputDeck1 Q-Chem input deck containing the first
    electronic state. 
    If value is a string, is treated as a filename for the QChemInput() output.
    If it is a QChemInput object, then it will be modified.
    If not specified, a new QChemInput object will be created.

    @param InputDeck2 Q-Chem input deck containing the second
    electronic state

    @param string NewFileName A new filename for the resultant QChemInput().
    (Default: tmp.in)

    @returns A Q-Chem input deck that is ready to be written and executed.
    
    NOTE: The file will NOT be written to disk. issue a .write() manually as
    necessary!   
    '''


    Q = QChemInput(NewFileName)
    del Q.jobs[0] #Delete empty job that is initialized

    if type(InputDeck1) == type('1'): #A string, treat as filename
        Q.ImportQChemOutput(InputDeck1)
    elif isinstance(InputDeck1, QChemInput):
        Q.jobs.append(deepcopy(InputDeck1.GetCurrentJob()))
    else:
        assert False, "I don't know what to do with InputDeck1 ="+\
            str(InputDeck1)

    Q.GetCurrentJob().append_input_section('rem', """\
symmetry                  off
sym_ignore                true
cdftci                    true
cdftci_stop               1
cdftci_skip_promolecules  true
print_input               false
qm_mm                     true
qmmm_print                true
skip_charge_self_interact true
print_orbitals            true
""")

    Q.GetCurrentJob().write_input_section('cdft', """\
0.
1. 0 0
0.
1. 0 0 s
---
0.
1. 0 0
0.
1. 0 0 s
""")

    if type(InputDeck2) == type('1'): #A string, treat as filename
        Q.ImportQChemOutput(InputDeck2)
    elif isinstance(InputDeck2, QChemInput):
        Q.jobs.append(deepcopy(InputDeck2.GetCurrentJob()))
    else:
        assert False, "I don't know what to do with InputDeck2 ="+\
str(InputDeck2)


    Q.GetCurrentJob().append_input_section('rem', """\
symmetry                  off
sym_ignore                true
cdftci                    true
cdftci_restart            1
cdftci_print              2
cdftci_skip_promolecules  true
print_input               false
qm_mm                     true
qmmm_print                true
skip_charge_self_interact true
""")

    Q.GetCurrentJob().write_input_section('cdft', """\
0.
1. 0 0
0.
1. 0 0 s
---
0.
1. 0 0
0.
1. 0 0 s
""")


    return Q


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


if __name__ == '__main__':
    import sys
    QCInput = QChemInputForElectrostaticEmbedding(*sys.argv[1:])
    if os.path.exists(QCInput.filename):
        print 'Output file', QCInput.filename, 'exists! Not overwriting.'
        print 'Writing to console instead.'
        print QCInput
    else:
        QCInput.write()
