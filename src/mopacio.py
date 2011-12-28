#!/usr/bin/env python
"""
Reads and writes MOPAC files

Contains routines and container classes that parse, represent and write
MOPAC input files, external parameter files, and output files.

Jiahao Chen <jiahao@mit.edu> 2011-04-20
"""

import logging, numpy, os, tempfile, unittest, UserDict
from subprocess import Popen, PIPE
import Stencil
from sysutils import logged_write
from copy import deepcopy

##Converts energies to gromacs units
eV = 96.485

##Convert forces to gromacs units
kcal_per_mol_ang = -41.8400

## List of valid MOPAC Parameter names
## as specified in @c src_subroutines/datin.F90:41
## in MOPAC 7.1 with some PM6 parameters appeneded
ValidMopacParameterNames = ['USS', 'UPP', 'UDD', 'ZS', 'ZP', 'ZD', 'BETAS',
    'BETAP', 'BETAD', 'GSS', 'GSP', 'GPP', 'GP2', 'HSP', 
    'AM1', 'EXPC', 'GAUSS', 'ALP', 'GSD', 'GPD', 'GDD',
    'ORBS', 'ZSN', 'ZPN', 'POCOR', 'BETASP', 'BETAPP','ZSN','ZDN','ZDP','G2SD', 'F0SD'] + \
     ['FN'+str(n)+str(m) for n in range(1,4) for m in range(1,11)]

## Maximum number of elements supported in MOPAC
## as specified in @c src_subroutines/datin.F90:46
## in MOPAC 7.1
MOPACMaxElements = 107

## List of known elements in MOPAC.
## as specified in @c src_subroutines/datin.F90:46
## in MOPAC 7.1
## @note In MOPAC's source code, elements 102--107 are special entities
## 102 CB is a capped bond, followed by four 'sparkles' ++ + -- -
## and TV for translation vectors
## @note Einsteinium ES is represented by XX
KnownElements = [ 'H', 'HE', 
                 'LI', 'BE',  'B',  'C',  'N',  'O',  'F', 'NE',
                 'NA', 'MG', 'AL', 'SI',  'P',  'S', 'CL', 'AR',
                  'K', 'CA',
                    'SC', 'TI',  'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN',
                             'GA', 'GE', 'AS', 'SE', 'BR', 'KR',
                 'RB', 'SR',
                     'Y', 'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD',
                             'IN', 'SN', 'SB', 'TE',  'I', 'XE',
                 'CS', 'BA',
                     'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU',
                     'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB',
                   'LU', 'HF', 'TA',  'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG',
                            'TL', 'PB', 'BI', 'PO', 'AT', 'RN',
                 'FR', 'RA',
                     'AC', 'TH', 'PA',  'U', 'NP', 'PU', 'AM',
                     'CM', 'BK', 'CF', 'XX', 'FM', 'MD', 'CB',
                            '++', '+', '--', '-' , 'TV'
                  ]


def vectorize(*args):
    """
    Generalized flattener.
    
    @param args,... Arbitrary number of iterables of any dimensionality
    @returns one-dimension numpy.ndarray()
    """
    buf = list()
    for arg in args:
        buf += list(numpy.ravel(numpy.asarray(arg)))
    return numpy.asarray(buf)


def isValidKey(key):
    """
    Checks that a given key is acceptable to MOPAC
    
    @param key 2-tuple containing key to validate. A valid key has a string
        containing a valid MOPAC keyword in index 0 and either an integer
        representing a valid atomic number in MOPAC or a string naming a
        valid element symbol in MOPAC.
    @returns @c True if key is a valid MOPAC keyword
    @throws @c KeyError if @c False
    """

    #The key must have two parts
    if len(key) != 2:
        raise KeyError(str(key)+' is not of length 2.')

    #The first part of the key must be a known MOPAC keyword.
    if key[0] not in ValidMopacParameterNames:
        if not (key[0].startswith("ALPB_") or key[0].startswith("XFAC_")):
            raise KeyError(str(key[0]) + ' is not a valid MOPAC keyword: ' 
                ', '.join(ValidMopacParameterNames))

    #The second part of the key must be an element: either the symbol
    #or a number
    try:
        isInteger = (repr(int(key[1])) == repr(key[1]))
        if isInteger:
            isValidInteger = 0 < int(key[1]) <= MOPACMaxElements
    except ValueError:
        isValidInteger = False

    if isValidInteger:
        return True
    
    try:
        isValidElementSymbol = key[1].upper() in KnownElements
    except (AttributeError, SyntaxError):
        isValidElementSymbol = False
    
    if not (isValidInteger or isValidElementSymbol):
        raise KeyError(str(key[1]) + ' lies outside valid range of [1, ' \
            + str(MOPACMaxElements) + '] or is not a valid element ' \
            + 'symbol: ' + ', '.join(KnownElements))
    
    return True



class MOPACParameters(UserDict.IterableUserDict):
    """
    A container for MOPAC parameters.
    """
    def __init__(self, params = None, filename = None):
        """
        @todo document
        @param params @todo document
        @param filename @todo document
        """
        UserDict.IterableUserDict.__init__(self)

        if params != None:
            for key, value in params.iteritems():
                assert isValidKey(key), 'Invalid key'
                self.data[key] = value

        ## Name of MOPAC parameters file
        self.filename = filename
        #If file exists, read it in
        try:
            self.read(self.filename)
        except (IOError, TypeError):
            pass


    def __iter__(self):
        """
        @todo document
        """
        for (paramname, element), value in self.data.iteritems():
            yield paramname, element, value


    def __delitem__(self, key):
        """
        @todo document
        @param key @todo document
        """
        if not isValidKey(key):
            raise KeyError, 'Invalid key: '+str(key)
        del self.data[key]

  
    def __getitem__(self, key):
        """
        @todo document
        @param key @todo document
        """
        if not isValidKey(key):
            raise KeyError, 'Invalid key: '+str(key)
        return self.data[key]
  

    def __setitem__(self, key, value):
        """
        @todo document
        @param key @todo document
        @param value @todo document
        """
        if not isValidKey(key):
            raise KeyError, 'Invalid key: '+str(key)
        self.data[key] = value


    def write(self, filename = None):
        """
        Writes parameters to a MOPAC external parameters file

        @param filename Name of external parameter file to write 
             or @c None, which will write parameter file to @c stdout
        """
        buf = []
        logging.info('Writing parameters')
        for paramname, element, value in self.__iter__():
            buf.append('   '.join(['', paramname, str(element), str(value)]))
            logging.debug(buf[-1])
        buf.append('')
        bufstr = '\n'.join(buf)
        
        #If a new filename is specified, overwrite
        if self.filename == None and filename != None:
            self.filename = filename

        if self.filename == None:
            print bufstr
            logging.info('Wrote MOPAC parameters to stdout')
        else:
            logged_write(self.filename, bufstr, logging)
            logging.info('Wrote MOPAC parameters to '+self.filename)


    def read(self, filename):
        """
        Reads parameters to a MOPAC external parameters file

        @param filename Name of external parameter file to read 
        """
        logging.info('Reading parameters from file '+filename)
        for line in open(filename):
            t = line.split()
            if len(t) == 0 or t[0].upper() == 'END':
                #Terminated by blank line or END
                break
            try:
                parname = t[0]
                element = t[1]
                if not isValidKey((parname, element)):
                    msg = 'Invalid key read:'+line
                    logging.critical(msg)
                    raise ValueError, msg
                value = float(t[2])
                logging.debug(line[:-1])
                self[parname, element] = value
            except (ValueError, IndexError), e:
                logging.critical('Invalid line read:')
                logging.critical(line[:-1])
                raise e
        logging.info('Finished reading parameters from '+filename)



class MOPACInput():
    """
    A container for MOPAC input decks
    """

    def __init__(self, filename = None, keywords = None):
        """
        Constructor.

        @param filename Name of MOPAC input file that it will represent.
               If this file exists, MOPACInput() will be initialized by parsing
               this file and keywords will be ignored.

        @param keywords Dictionary of keywords. Default: None
        """

        ## Name of MOPAC input file represented
        self.filename = filename

        ## First comment line
        self.comment1 = ''

        ## Second comment line
        self.comment2 = ''
        
        if keywords == None:
            self.keywords = dict()
        else:
            ## Dictionary of keywords
            self.keywords = keywords

        ## List of geometries
        ## @todo document
        self.geometry = []

        ## List of geometry flags of optimizable coordinates
        ## @todo document
        self.geomoptflags = []

        try:
            self.read(filename)
        except (TypeError, IOError):
            pass



    def read(self, filename):
        """
        Parses MOPAC input deck for keywords, comments and geometry

        @param filename Name of MOPAC input deck to read

        @todo Currently only supports input decks with only one line of
              keywords. Want to support continued input lines.
        @todo Support specification of other data
        """
        state = 'commands'
        for line in open(filename):
            t = line.split()

            if state == 'commands':
                for keywordvalpair in t:
                    logging.debug('Found keyword: '+keywordvalpair)
                    try:
                        keyword, value = keywordvalpair.split('=')
                    except ValueError:
                        keyword, value = keywordvalpair, None
                    self.keywords[keyword] = value
                state = 'comment1'

            elif state == 'comment1':
                self.comment1 = line[:-1]
                logging.debug('Comment line 1: '+self.comment1)
                state = 'comment2'

            elif state == 'comment2':
                self.comment2 = line[:-1]
                logging.debug('Comment line 2: '+self.comment2)
                logging.debug('Reading geometry:')
                state = 'geometry'

            elif state == 'geometry':
                try:
                    element = t[0]
                    x, y, z = float(t[1]), float(t[3]), float(t[5])
                    flagx, flagy, flagz = int(t[2]), int(t[4]), int(t[6])
                    logging.debug(line[:-1])
                    self.geometry.append((element, x, y, z))
                    self.geomoptflags.append((flagx, flagy, flagz))
                except (IndexError, ValueError), e:
                    logging.critical('ONLY CARTESIAN GEOMETRIES SUPPORTED')
                    raise e



    def write(self, filename = None):
        """
        Writes MOPAC input deck

        @param filename Name of input deck to write
        """

        #Assemble keywords
        logging.debug('Assembling MOPAC input')
        kwbuf = []
        for key, val in self.keywords.items():
            if val == None:
                kwbuf.append(key)
            else:
                kwbuf.append(key+'='+str(val))
            logging.debug('MOPAC keyword '+kwbuf[-1])
        keywords = ' '.join(kwbuf)

        #Now assemble comments and geometry
        buf = [keywords, self.comment1, self.comment2]
        logging.debug('Comment line 1'+self.comment1)
        logging.debug('Comment line 2'+self.comment2)

        logging.debug('Cartesian geometry')
        for i, line in enumerate(self.geometry):
            try:
                flagline = self.geomoptflags[i]
            except IndexError:
                # The default is to assume all 1s, i.e. all coordinates are
                # optimizable and hence will have derivatives calculated and
                # printed
                flagline = (1, )*3

            geometryline = "%s %f %d %f %d %f %d" % (line[0], line[1], \
                    flagline[0], line[2], flagline[1], line[3], flagline[2])
            buf.append(geometryline)
            logging.debug(geometryline)

        bufstr = '\n'.join(buf)

        #If a new filename is specified, overwrite
        if self.filename == None and filename != None:
            self.filename = filename

        if self.filename == None:
            print bufstr
            logging.info('Wrote MOPAC input to stdout')
        else:
            f = open(self.filename, 'w')
            f.write(bufstr)
            f.close()
            logging.info('Wrote MOPAC input to '+self.filename)



    def getfilenameroot(self, filename = None):
        """Get the root part of the filename
        
        MOPAC 2009's logic is:

        If filename ends in .dat or .mop (in any case), root is the
        part without the extension, otherwise the entire input filename
        is the root.
        """
        if filename is None:
            filename = self.filename

        if filename[-4:].lower() == '.dat':
            return filename[:-4]
        elif filename[-4:].lower() == '.mop':
            return filename[:-4]
        else:
            return filename

    def execute(self, dowrite = True, cmd = 'mopac'):
        """
        Executes the current input deck.
        
        @param cmd Name of MOPAC executable. Default: 'mopac'
        @param dowrite if @c True, also writes input deck using write().
               Default: True
        @returns A MOPACOutput() containing the parsed output from MOPAC
        """
        if dowrite:
            self.write()
        cmdstring = ' '.join([cmd, self.filename])

        p = Popen(cmdstring, shell=True, stderr=PIPE)

        p.wait()

        #Looks for errors
        #By default MOPAC2009 prints something like
        #          Mon Jun 20 20:22:32 2011  Job: 'si' started successfully
        #to stderr
        error_msg = p.stderr.readlines()
        assert not error_msg or 'started successfully' in error_msg[0], \
                '\n'.join(error_msg)

        return MOPACOutput(self.getfilenameroot()+'.out')



class MOPACPeriodicInput(MOPACInput):
    """
    MOPAC Input container for periodic systems

    Identical to MOPACInput() except that for convenience, additional keywords
    are turned on during initialization.

    The extra things that must be specified are

    Translation vectors: specify in geometry record as TV 
    DEBUG DCART LARGE: for gradient on central unit cell to be printed

    @see http://openmopac.net/manual/Solids_derivatives.html
    """

    def __init__(self, filename = None, keywords = None):
        """
        Constructor.

        @param filename Name of MOPAC input file that it will represent.
               If this file exists, MOPACPeriodicInput() will be initialized
               by parsing this file and keywords will be ignored.
        @param keywords Dictionary of keywords. Default: None
               
        @note The keywords DEBUG DCART LARGE will *always* be added.
        """

        MOPACInput.__init__(self, filename, keywords)
        self.keywords['DEBUG'] = None 
        self.keywords['DCART'] = None 
        self.keywords['LARGE'] = None 



class _MOPACMatrixData:
    """Data container to aid parsing of matrix output in MOPACOutput()"""
    def __init__(self):

        ## A dictionary mapping header labels to internal indices
        self.Labels = {}

        ## This will be a dictionary mapping (i,j) indices to a value
        self.MatrixElements = {}

        ## Keeps track of the current column labels
        self.ColLabels = []

    def get_asarray(self):
        """Returns self in numpy array format"""
        N = len(self.Labels)
        M = numpy.empty((N, N))
        for (i, j), x in self.MatrixElements.items():
            M[i, j] = M[j, i] = x

        return M

    def get_labels(self):
        """
        Returns labels in numpy array format
        """

        N = len(self.Labels)
        datatype = numpy.dtype([('AO', 'S2'), ('Element', 'S2'), 
                                ('Index', numpy.int16)])

        Labels = numpy.empty((N, ), dtype = datatype)
        for label, i in self.Labels.items():
            Labels[i] = label

        return Labels


class _MOPACOutputParser:
    """Parser for MOPAC Output file

    This class encapsulates the finite state machine used for the parsing
    of the output file.
    """
    def __init__(self, filename = None, version = '2009'):
        """Initializer
        
        @param filename Name of MOPAC Output file to parse. If not @c None.
        launches parse() immediately to being parsing. Default: @c None.
        
        @param version MOPAC Version to parse. Valid values: '7.1' or
        '2009'. Default: '2009'.
        """
        ## Internal state of the parser 
        self.state = 'seek'
        ## Dictionary of all data found by the parser 
        self.data = dict()
        ## Current data set
        self.thisdata = None
        ## Name of current data set
        self.thisdataname = None
        ## Number of lines to skip
        self.skiplines = 0
        ## MOPAC Version. This is used when there are version-specific
        ##things to worry about during parsing. Valid values: '7.1' or '2009'
        self.MOPAC_version = version

        if filename is not None:
            logging.info('Parsing %s', filename)
            self.parse(filename)

    def parse(self, filename):
        """Start parser
        
        @param filename Name of MOPAC output file to parse
        
        @note parse() will farm out calls to internal private functions
        based on the current state. For example, if self.state == 'seek',
        parse() will call _seek(line) using the current @c line."""
        logging.debug('Opening file %s', filename)
        
        for line in open(filename):
            #print self.state, line,
            if self.skiplines > 0:
                self.skiplines -= 1
            else:
                logging.debug('%s: %s', self.state, line)
                getattr(self, '_' + self.state)(line)
        
        if 'NormalTermination' not in self.data:
            logging.warning('MOPAC calculation did not terminate normally')

    ### The following private methods are designed to be called by the
    ### parser only!
    def _seek(self, line):
        """
        Looks for recognized titles that caption the beginning of a
        recognized data set.
        """
        if 'MATRIX' in line and line[:2] != ' *':
            logging.debug('Found matrix output')
            self.thisdataname = line.strip()
            self.thisdata = _MOPACMatrixData()
            self.skiplines = 2
            self.state = 'read_matrix'

        elif '== MOPAC DONE ==' in line:
            self.data['NormalTermination'] = True
            logging.debug('MOPAC calculation terminated normally')

        elif 'CARTESIAN COORDINATES' in line:
            logging.debug('Found Cartesian coordinate block')
            self.state = 'geometry'
            self.thisdata = []
            self.skiplines = 3

        elif 'PARAMETER TYPE      ELEMENT    PARAMETER' in line \
              or 'Parameter Type  Element    Parameter' in line:  #MOPAC 2009
            logging.debug('Found parameter specification block')
            self.state = 'parameters'

        elif 'CARTESIAN COORDINATE DERIVATIVES' in line:
            logging.debug('Found cartesian derivatives block')
            self.state = 'cartesian_derivative'
            self.thisdata = []
                    
            if self.MOPAC_version == '7.1':
                #For some reason, MOPAC prints out the derivatives
                #twice, so we will skip the first of the repeated blocks
                self.skiplines = 8 + self.data['NumAtoms']
            elif self.MOPAC_version == '2009':
                self.skiplines = 3

        elif 'Central Unit Cell Derivatives' in line:
            logging.debug('Found central unit cell derivatives')
            self.state = 'cuc_derivative'
            self.thisdata = []

        elif 'TOTAL ENERGY' in line:
            t = line.split()
            thisenergy = float(t[-2])
            thiseunit = t[-1]
            if thiseunit == 'EV':
                thisenergy *= eV
            else:
                logging.error('Unknown energy unit: '+t[-1])
                logging.error(line[-1])
                raise ValueError

            if 'Energies' not in self.data:
                self.data['Energies'] = list()

            self.data['Energies'].append(thisenergy)
            logging.info('Read in total energy = %f kJ/mol', thisenergy)

    def _read_matrix(self, line):
        """
        Parse in a line of matrix output
        MOPAC outputs matrices in lower triangular format
        """

        t = line.split()
        try:
            float(t[3])
            LineType = 'data'
        except (ValueError, IndexError):
            if len(t) >= 2:
                LineType = 'header'
            else: #Ignore
                return

        if LineType == 'header':
            # I really hate this but it has to be done this way
            # MOPAC is dumb and lets the numbers run on
            # So one has ' S Si  1' but also 'PX Tc242'
            # therefore instead of tokenizing I've resorted to
            # manual character counts *sigh*
            # To make things worse the row labels are 'PZ Cl   2'
            # which is INCONSISTENT /rant

            length_headertoken = 11
            startpos = 12
            numheaders = (len(line)-startpos)/(1.0 * length_headertoken)

            assert numheaders == int(numheaders), 'Format of header not \
recognized: '+line

            self.thisdata.ColLabels = list()

            for i in range(int(numheaders)):
                thisheader = line[(14+i*length_headertoken):\
                                  (22+i*length_headertoken)]
                thisheader = (thisheader[:3], thisheader[3:5], thisheader[5:])
                thisheader = (thisheader[0].strip(), thisheader[1].strip(), 
                              int(thisheader[2]))
                #This should correctly extract 'PX Tc242' --> ('PX', 'Tc', 242)

                try:
                    colidx = self.thisdata.Labels[thisheader]
                except KeyError: #Not already in there
                    colidx = len(self.thisdata.Labels)
                    self.thisdata.Labels[thisheader] = colidx

                self.thisdata.ColLabels.append((thisheader, colidx))

            #print 'New column labels'
            #print self.ColLabels

        elif LineType == 'data':
            try:
                thisheader = (t[0], t[1], int(t[2]))
            except (ValueError, IndexError):
                #Assume done
                self.data[self.thisdataname] = self.thisdata.get_asarray()
                self.data[self.thisdataname+' LABELS'] = \
                        self.thisdata.get_labels()
                self.state = 'seek'
                return

            try:
                rowidx = self.thisdata.Labels[thisheader]
            except KeyError: #Not already in there
                rowidx = len(self.thisdata.Labels)
                self.thisdata.Labels[thisheader] = rowidx

            for i, x in enumerate(t[3:]):
                colidx = self.thisdata.ColLabels[i][1]
                self.thisdata.MatrixElements[rowidx, colidx] = float(x)

            #print rowidx, [float(x) for x in t[3:]]

    def _geometry(self, line):
        """Reading in a geometry line"""
        # Blank line terminates geometry in MOPAC 7.1
        # Parameters teminates geometry in MOPAC 2009
        t = line.split()
        if len(t) == 0 or t[0] == 'Parameters':

            if 'Geometries' not in self.data:
                self.data['Geometries'] = []

            self.data['Geometries'].append(self.thisdata)
            assert self.data['NumAtoms'] == len(self.thisdata), \
                    'Did not read in the correct number of atoms'

            logging.info('Read in %s atoms', self.data['NumAtoms'])
            self.state = 'seek'
            return

        try:
            self.data['NumAtoms'] = int(t[0])
            if len(t) == 5:
                atom = t[1]
                x, y, z = float(t[2]), float(t[3]), float(t[4])
            else:
                templine = line.replace("-"," -")
                print templine,
                t = line.split()
                atom = t[1]
                x, y, z = float(t[2]), float(t[3]), float(t[4])
            self.thisdata.append((atom, x, y, z))
            logging.debug(line[:-1])
        except (IndexError, ValueError), e:
            logging.error('Invalid line encountered in geometry:')
            logging.error(line[:-1])
            raise e

    def _parameters(self, line):
        """Reading in parameters"""
        t = line.split()
        if len(t) == 0: #Blank line terminates block
            self.state = 'seek'
            logging.debug('Finished reading in parameters')
            return

        #Initialize
        if 'Parameters' not in self.data:
            self.data['Parameters'] = dict()

        try:
            paramname = t[0]
            element = t[1]
            value = float(t[2])
            self.data['Parameters'][paramname, element] = value
            logging.debug(line[:-1])
        except (IndexError, ValueError), e:
            logging.error('Invalid line encountered in geometry:')
            logging.error(line[:-1])
            raise e

    def _cartesian_derivative(self, line):
        """Read in a line of cartesian derivatives"""
        #There is a workaround MOPAC2009 insanity where for a periodic
        #system, it prints an insanely long CARTESIAN COORDINATE
        #DERIVATIVES block
        #
        #To avoid this, check length of derivative.
        #If too long, don't store

        if len(self.thisdata) > self.data['NumAtoms']:
            self.state = 'seek'
            return

        #Convert force units to gromacs units
        unit = kcal_per_mol_ang
        
        t = line.split()
        if len(t) == 5 or len(t) == 4:
            #For example 
            # 3     H   854.785342    23.681706  1895.110923
            try:
                x, y, z = float(t[-3]), float(t[-2]), float(t[-1])
            except ValueError, e:
                logging.error('Error parsing geometry: '+line)
                raise e
            self.thisdata.append((x*unit, y*unit, z*unit))
            
        elif len(t) == 7:
            #For example
            #39Si     0.000062     0.000294     0.000197  -2   0   2
            #This format is only produced by when there are
            #periodic boundary conditions in MOPAC2009. Ignore.
            self.state = 'seek'

        else:
            #Assume we are done
            if 'Cartesian derivatives' not in self.data:
                self.data['Cartesian derivatives'] = list()

            self.data['Cartesian derivatives'].append(self.thisdata)
            logging.debug('Finished reading in derivative')
            self.state = 'seek'
        
    def _cuc_derivative(self, line):
        """Cartesian derivative for central unit cell 
        for use with periodic boundary conditions"""

        #Check for end of block
        if 'Fractional Unit Cell Derivatives' in line:
            if 'Cartesian derivatives' not in self.data:
                self.data['Cartesian derivatives'] = list()
            self.data['Cartesian derivatives'].append(self.thisdata)
            logging.debug('Finished reading in derivative')
            self.state = 'seek'
            return

        #Convert force units to gromacs units
        unit = kcal_per_mol_ang
        #print line
        t = line.split()
        if len(t) == 5 or len(t) == 4:
            try:
                x, y, z = float(t[-3]), float(t[-2]), float(t[-1])
            except ValueError, e:
                logging.error('Error parsing geometry: %s', line)
                raise e
            self.thisdata.append((x*unit, y*unit, z*unit))
        else:
            #Assume we are done
            if 'Cartesian derivatives' not in self.data:
                self.data['Cartesian derivatives'] = list()
            self.data['Cartesian derivatives'].append(self.thisdata)
            logging.debug('Finished reading in derivative')
            self.state = 'seek'



class MOPACOutput():
    """A container for parsed MOPAC output files"""
    def __init__(self, filename = None):
        """
        Constructor

        @param filename Name of MOPAC output file represented.
               If the file exists, uses read() to initialize data.
               Default: None
        """

        ## All the parsed data, ever.
        self.data = dict()
        ## Name of output file
        self.filename = filename
        #If file exists, read it in
        self.read(filename)

    def read(self, filename):
        """
        Reads and parses MOPAC output file

        @param filename Name of file to parse.
        """
        logging.info('Calling output parser')  
        self.data = _MOPACOutputParser(filename).data

    def get_data(self, dataname = None):
        """Returns data found by parser.

        @param dataname Name of data set to return, or @c None.
            If @c None, returns all data. Default: @c None.

        @returns Data set with name @c dataname. If not found, returns
            @c None. If @c dataname is @c None, returns the entire
            collection of data in self.data.
        """

        if dataname is None:
            return self.data

        if dataname not in self.data:
            logging.error('Did not parse requested data: '+dataname)
            return None
        else:
            return self.data[dataname]

    def CalculateSpinPopulations(self, doPrint = True):
        """
        Calculates atomic spin populations using Coulson-type
        population analysis.

        For an atom @f$ i @f$, the Coulson spin population is defined as
        @f[
            S_i = @sum_{@mu @in i} P^s_{@mu@mu}

        @f]

        where @f$ @mu @f$ indexes atomic orbitals and @f$ P^s @f$ is the
        spin density matrix.

        @param doPrint Boolean specifying whether or not to print the
        calculated atomic populations to stdout. Default: @c True
        """

        assert 'SPIN DENSITY MATRIX' in self.data, \
            'Insufficient data to calculate spin populations'

        assert 'SPIN DENSITY MATRIX LABELS' in self.data, \
            'Insufficient data2 to calculate spin populations'

        populations = numpy.diag(self.data['SPIN DENSITY MATRIX'])
        labels = self.data['SPIN DENSITY MATRIX LABELS']

        numAtoms = labels[-1][2] #Should be max over labels[:][2] but we assume
                                 #here output in ascending order

        atomicPopulations = numpy.zeros((numAtoms, ))

        for i, population in enumerate(populations):
            atomIndex = labels[i][2] - 1
            atomicPopulations[atomIndex] += population

        if doPrint:
            print 'Calculation of Coulson spin populations requested'
            print
            print 'Atom ID    Population'
            print '---------------------'

        for i, p in enumerate(atomicPopulations):
            atomSymbol = [label[1] for label in labels if label[2] == i+1][0]
            if doPrint:
                print '%2s %4d     %9.6f' % (atomSymbol, i+1, p)

        return atomicPopulations



def MOPACGetEnergyAndForce(paramvals, paramkeys, geometry, keywords,
                           parfilename, inputfilename):
    """
    Calculates energy and force using MOPAC

    @param paramvals A one-dimensional iterable of parameter values
    @param paramkeys A two-dimensional iterable of MOPAC keyword semantics
           of the form [(ParameterName1, ElementSymbolOrNumber1), 1<-->2, ...]
           over parameters

    @param geometry  A two-dimensional iterable containing a Cartesian geometry
           of the form [(ElementSymbolOrNumber1, x1, y1, z1), 1<-->2, ...]
           over atoms

    @param keywords  A dictionary of keywords to be specified in the MOPAC
           input deck whose values are the arguments to be passed to the
           keywords.

           MOPAC keywords themselves can have zero or one arguments.
           For example, ('PM3', None) becomes the input deck keyword PM3
                    and ('GNORM', 1e-3) becomes GNORM=0.001

    @param parfilename Name of the parameter file to be generated
    @param inputfilename Name of the MOPAC input deck to be generated
    """

    #Assemble internal representation of parameters
    parameters = dict()
    for pkey, pval in zip(paramkeys, paramvals):
        parameters[pkey] = pval

    PM3 = MOPACParameters(parameters)
    PM3.write(parfilename)
    f = MOPACInput(inputfilename)
    f.geometry = geometry
    f.keywords = keywords
    f.write()
    output = f.execute(cmd = 'mopac2009')

    #To avoid retarded error arising from repeated execution,
    #clean up archive file
    
    arc_filename = f.getfilenameroot() + '.arc'
    if os.path.exists(arc_filename):
        os.unlink(arc_filename)

    #Take the first energy and the first cartesian derivative block
    energy = output.get_data('Energies')[0]
    force = output.get_data('Cartesian derivatives')[0]

    #Vectorize the output for convenience in calling from
    #MOPACGetEnergyForceAndDerivs()
    return vectorize(energy, force)



def MOPACGetEnergyForceAndDerivs(paramvals, paramkeys, geometry, keywords,
     parfilename, inputfilename, Stc = None, step = 0.001):
    """
    Calculates energy and force using MOPACGetEnergyAndForce(), and also
    calculates their derivatives with respect parameters using a numerical
    finite difference stencil

    @param paramvals A one-dimensional iterable of parameter values
    @param paramkeys A two-dimensional iterable of MOPAC keyword semantics
           of the form [(ParameterName1, ElementSymbolOrNumber1), 1<-->2, ...]
           over parameters

    @param geometry  A two-dimensional iterable containing a Cartesian geometry
           of the form [(ElementSymbolOrNumber1, x1, y1, z1), 1<-->2, ...]
           over atoms

    @param keywords  A dictionary of keywords to be specified in the MOPAC
           input deck whose values are the arguments to be passed to the
           keywords.

           MOPAC keywords themselves can have zero or one arguments.
           For example, ('PM3', None) becomes the input deck keyword PM3
                    and ('GNORM', 1e-3) becomes GNORM=0.001

    @param parfilename Name of the parameter file to be generated
    @param inputfilename Name of the MOPAC input deck to be generated
    @param Stc a Stencil.Stencil to be used. (Default:
           Stencil.FirstOrderCentralDifferenceStencil )
    @param step Magnitude of finite difference step to be used in Stc.
           (Default: 0.001)
    """


    ### MGW: Added density saving/loading to hopefully speed up finite difference
    keys_denout = deepcopy(keywords)
    keys_denout['DENOUT'] = None
    data = MOPACGetEnergyAndForce(paramvals, paramkeys, geometry, keys_denout,
                                  parfilename, inputfilename)
    Energy, Force = data[0], data[1:]

    #Assign default stencil
    if Stc == None:
        Stc = Stencil.FirstOrderCentralDifferenceStencil()
    
    keys_oldens = deepcopy(keywords)
    keys_oldens['OLDENS'] = None
    derivdata = Stc.ApplyToFunction(MOPACGetEnergyAndForce,
                   paramvals, step, paramkeys, geometry, keys_oldens, parfilename,
                   inputfilename)

    EnergyDeriv = derivdata[0, :]
    ForceDeriv  = derivdata[1:, :]

    return Energy, Force, EnergyDeriv, ForceDeriv



class TestMOPAC(unittest.TestCase):
    """
    Unit test class for MOPAC 2009
    """

    def test_methane(self):
        """
        A simple unit test using methane.
        """
        
        parfilename = tempfile.mktemp(suffix='.par')
        inputfilename = tempfile.mktemp(suffix='.mop')
    
        #INPUTS
        geometry = [
            ('C', -1.10303,  0.80250,  0.0),
            ('H', -0.03303,  0.80250,  0.0),
            ('H', -1.54753,  0.78560, -0.97316), 
            ('H', -1.45795,  1.65691,  0.53751), 
            ('H', -1.39785, -0.11001,  0.47467)
        ]
    
        #parameter values: an ordered vector of floats
        paramvals = [0.967807, 2.707807, 0.050107, 6.002979, -1.060329,
                     1.537465, 2.29098, -47.27032, -9.802755, 1.565085,
                    -5.626512, 10.796292, 0.892488, 1.570189, 3.356386,
                     9.042566, 1.12875, 1.642214, 0.050733, 10.265027,
                     6.003788, 14.794208, -11.910015, -36.266918, -13.073321,
                     1.842345, 11.200708, 6.003165, 5.096282]
    
        #A static list of semantic qualifies for the vector of parameters
        paramkeys = [('ZS', 'H'), ('ALP', 'C'), ('FN11', 'C'), ('FN22', 'C'),
                     ('FN21', 'H'), ('FN13', 'H'), ('HSP', 'C'), ('USS', 'C'),
                     ('BETAP', 'C'), ('ZP', 'C'), ('BETAS', 'H'), ('GPP', 'C'),
                     ('FN23', 'C'), ('FN23', 'H'), ('ALP', 'H'), ('GP2', 'C'),
                     ('FN11', 'H'), ('FN13', 'C'), ('FN21', 'C'), ('GSP', 'C'),
                     ('FN22', 'H'), ('GSS', 'H'), ('BETAS', 'C'), ('UPP', 'C'),
                     ('USS', 'H'), ('ZS', 'C'), ('GSS', 'C'), ('FN12', 'C'),
                     ('FN12', 'H')]
    
        #Default list of keywords
        defaultkeywords = {
            'GNORM': '9E99',         ### Gradient norm threshold for auto-
                                     ### triggering geometry optimization
            'PM3': None,             ### Use PM3 Hamiltonian
            'GEO-OK': None,          ### Turn off geometry sanity checks
            'DCART': None,           ### Print out Cartesian energy gradients
            'EXTERNAL': parfilename, ### Name of external parameter file
        #    'ANALYT': None          ### Use analytic gradients
            }
    
        #########
    
        #Now do the actual calculations
        Energy, Force, EnergyDeriv, ForceDeriv = MOPACGetEnergyForceAndDerivs(\
            paramvals, paramkeys, geometry, defaultkeywords, parfilename,
            inputfilename)
    
        #Clean up
        inputroot = MOPACInput().getfilenameroot(inputfilename)
        for filename in (inputfilename, parfilename, inputroot+'.out',
             inputroot+'.temp', inputroot+'.log', inputroot+'.end',
             inputroot+'.brz', inputroot+'.arc'):
            if os.path.exists(filename): 
                os.unlink(filename)
    
        logging.info('Energy = '+ str(Energy))
        logging.info('Force  = '+ str(Force))
        logging.info('Energy derivative = '+ str(EnergyDeriv))
        logging.info('Force  derivative = '+ str( ForceDeriv))
    
        #print repr(Energy), repr(Force), repr(EnergyDeriv), repr(ForceDeriv)
        
        savedEnergy = -2103.1829245500003

        savedForce = numpy.array([\
         1165.98615792,   3626.16393232,  -1880.02408336,  86990.61615   ,
          379.30311408,   -414.5461176 , -35763.40864504,   -990.83973392,
       -79289.695328  , -28811.76674368,  70305.65099208,  43614.03633424,
       -23581.4269192 , -73320.27830456,  37970.22919472])

        savedEnergyDeriv = numpy.array([\
        -589.040925,  -843.2789  ,  1418.811925,  -226.257325,
        -483.872275,   244.589475,  -205.030625,   105.651075,
         159.682675,  1300.6178  ,   252.308275,   -24.603675,
        -976.4282  ,  -446.72555 ,  -463.128   ,   364.230875,
        2513.916675,   519.571725,   -23.638825,   307.304725,
        -496.415325,    58.373425,   103.721375,   280.288925,
         385.457575,   837.972225,    13.025475,  1194.001875,   869.812275]) 
        savedForceDeriv = numpy.array([[\
         -2.49010760e+02,  -1.05750600e+02,   7.03958000e+01,
          1.04851040e+02,  -1.57423000e+02,   5.80739200e+01,
          0.00000000e+00,  -5.14632000e+00,   7.59396000e+00,
         -5.33041600e+01,   8.97468000e+00,  -2.30120000e+00,
         -5.39526800e+01,  -5.77810400e+01,  -6.82828800e+01,
         -5.96220000e+00,   1.98196080e+02,   8.35963200e+01,
         -5.77392000e+00,  -1.76355600e+01,   2.47902000e+01,
          2.46856000e+00,   8.28432000e+00,  -3.24260000e+00,
          8.66088000e+00,   3.52920400e+01,   6.90360000e-01,
          2.07317200e+02,  -9.22990400e+01],
       [ -7.81926840e+02,  -3.27377080e+02,   2.17944560e+02,
          3.24594720e+02,  -4.87707960e+02,   1.79870160e+02,
          1.96648000e+00,  -1.43302000e+01,   2.40370800e+01,
         -1.52025640e+02,   2.46228400e+01,  -7.57304000e+00,
         -1.66899760e+02,  -1.78845080e+02,  -2.11417520e+02,
         -2.60244800e+01,   6.13792800e+02,   2.58822240e+02,
         -1.78447600e+01,  -5.47476400e+01,   7.69856000e+01,
          9.22572000e+00,   2.40998400e+01,  -1.37235200e+01,
          2.85767200e+01,   1.01964080e+02,   3.17984000e+00,
          6.41511800e+02,  -2.85516160e+02],
       [  4.01371120e+02,   1.70435240e+02,  -1.13470080e+02,
         -1.68949920e+02,   2.54094320e+02,  -9.36797600e+01,
         -4.18400001e-02,   8.28432000e+00,  -1.22382000e+01,
          8.57092400e+01,  -1.43929600e+01,   3.72376000e+00,
          8.68180000e+01,   9.30730800e+01,   1.10081040e+02,
          9.72780000e+00,  -3.19657600e+02,  -1.34766640e+02,
          9.28848000e+00,   2.84093600e+01,  -4.02082400e+01,
         -4.01664000e+00,  -1.33051200e+01,   5.29276000e+00,
         -1.39745600e+01,  -5.66304400e+01,  -1.12968000e+00,
         -3.33799520e+02,   1.48490160e+02],
       [ -7.01309528e+03,  -4.67108036e+03,   3.10946512e+03,
          4.63120684e+03,  -5.59082816e+03,   2.70121132e+03,
          2.85976400e+01,  -1.46314480e+02,   1.42611640e+02,
         -3.86231316e+03,   4.99276720e+02,   1.56900000e+00,
         -2.38082152e+03,  -3.06300180e+03,  -3.55627448e+03,
         -7.55839600e+01,   9.08474012e+03,   3.69298668e+03,
         -2.54449960e+02,  -4.58650080e+02,   4.03461028e+03,
          1.22486600e+02,   2.81625040e+02,  -3.06478000e+01,
          1.76836760e+02,  -3.44740680e+02,   9.49558800e+01,
          9.15226988e+03,   1.92440988e+03],
       [  1.17256600e+02,   0.00000000e+00,   0.00000000e+00,
          0.00000000e+00,  -2.05225200e+01,   9.68596000e+00,
          2.71960000e-01,   2.74052000e+00,  -5.96220000e+00,
         -6.40570400e+01,   2.61500000e+00,   2.32212000e+00,
          0.00000000e+00,  -2.98319200e+01,  -2.65056400e+01,
          8.05420000e+00,   6.40779600e+01,   0.00000000e+00,
          0.00000000e+00,   1.14432400e+01,   3.67355200e+01,
         -3.13800000e-01,  -4.12124000e+00,   4.20492000e+00,
         -7.17556000e+00,  -3.57732000e+01,   2.15476000e+00,
          0.00000000e+00,   2.77608400e+01],
       [ -2.13279400e+02,   0.00000000e+00,   0.00000000e+00,
          0.00000000e+00,   8.14624800e+01,  -1.40164000e+01,
          4.18400000e-02,   8.36800000e-02,   1.47276800e+01,
          1.52025640e+02,   1.46440000e+00,  -6.52704000e+00,
          0.00000000e+00,   4.31370400e+01,   3.66727600e+01,
         -2.55851600e+01,  -1.19076640e+02,   0.00000000e+00,
          0.00000000e+00,  -1.38908800e+01,  -1.09181480e+02,
          1.12968000e+00,   2.05016000e+00,  -1.34306400e+01,
          1.34306400e+01,   2.46437600e+01,  -4.18400000e-01,
          0.00000000e+00,   1.12549600e+01],
       [  3.10521836e+03,   1.94020448e+03,  -1.29155896e+03,
         -1.92371952e+03,   2.24787492e+03,  -1.10930392e+03,
         -1.19244000e+01,   6.06680000e+01,  -7.25505600e+01,
          1.46705684e+03,  -2.08760680e+02,   5.27184000e+00,
          9.89076680e+02,   1.23354780e+03,   1.44419128e+03,
          5.44966000e+01,  -3.66585344e+03,  -1.53389624e+03,
          1.05687840e+02,   2.02944920e+02,  -1.57699144e+03,
         -5.18816000e+01,  -1.18783760e+02,   2.48738800e+01,
         -8.55628000e+01,   1.21273240e+02,  -3.90785600e+01,
         -3.80195896e+03,  -8.10064240e+02],
       [  2.30120000e+02,   7.37639200e+01,  -4.90992400e+01,
         -7.31363200e+01,   6.69021600e+01,  -3.28234800e+01,
         -1.88280000e-01,   5.04172000e+00,  -8.34708000e+00,
         -4.60240000e+00,  -5.29276000e+00,   2.36396000e+00,
          3.76141600e+01,   1.80958000e+01,   2.92670800e+01,
          9.49768000e+00,  -7.81780400e+01,  -5.83040400e+01,
          4.03756000e+00,   1.88489200e+01,  -2.59826400e+01,
         -2.23844000e+00,  -8.59812000e+00,   4.83252000e+00,
         -1.01043600e+01,  -3.05850400e+01,   6.48520000e-01,
         -1.44536280e+02,  -2.82420000e+00],
       [  6.28825912e+03,   4.24803612e+03,  -2.82784008e+03,
         -4.21186544e+03,   5.11791064e+03,  -2.46238860e+03,
         -2.59826400e+01,   1.33072120e+02,  -1.23532600e+02,
          3.57639952e+03,  -4.53503760e+02,  -4.14216000e+00,
          2.16538736e+03,   2.80390760e+03,   3.24969188e+03,
          5.79693200e+01,  -8.31176704e+03,  -3.35845496e+03,
          2.31396120e+02,   4.11224440e+02,  -3.71472256e+03,
         -1.10917840e+02,  -2.55224000e+02,   2.22379600e+01,
         -1.55142720e+02,   3.24008960e+02,  -8.65251200e+01,
         -8.32385880e+03,  -1.74569032e+03],
       [  2.30942156e+03,   1.54966992e+03,  -1.03158612e+03,
         -1.53634388e+03,   1.87754908e+03,  -8.95878080e+02,
         -9.16296000e+00,   5.06264000e+01,  -4.47060400e+01,
          1.30743724e+03,  -1.63803600e+02,  -1.86188000e+00,
          7.89688160e+02,   1.01503840e+03,   1.17823532e+03,
          1.83468400e+01,  -3.02383956e+03,  -1.22524256e+03,
          8.44122000e+01,   1.53720160e+02,  -1.36791696e+03,
         -4.01873200e+01,  -9.54998000e+01,   6.58980000e+00,
         -5.74672400e+01,   1.00834400e+02,  -3.01038800e+01,
         -3.03588948e+03,  -6.20131560e+02],
       [ -5.62105756e+03,  -3.73001508e+03,   2.48299480e+03,
          3.69813300e+03,  -4.42775984e+03,   2.15354664e+03,
          2.29910800e+01,  -1.14223200e+02,   1.17695920e+02,
         -3.05030336e+03,   4.01915040e+02,  -8.36799998e-01,
         -1.90110500e+03,  -2.43897912e+03,  -2.83597796e+03,
         -6.86385200e+01,   7.22773448e+03,   2.94896688e+03,
         -2.03195960e+02,  -3.63987080e+02,   3.34987776e+03,
          9.80729600e+01,   2.22212240e+02,  -2.88277600e+01,
          1.42925440e+02,  -2.96122600e+02,   7.75713600e+01,
          7.30823464e+03,   1.72083736e+03],
       [ -3.49619224e+03,  -2.34665916e+03,   1.56211732e+03,
          2.32655504e+03,  -2.84407400e+03,   1.35653648e+03,
          1.38699600e+01,  -7.68182400e+01,   6.75297600e+01,
         -1.98158424e+03,   2.47943840e+02,   2.92880000e+00,
         -1.19597548e+03,  -1.53720160e+03,  -1.78412036e+03,
         -2.73006000e+01,   4.57942984e+03,   1.85531112e+03,
         -1.27821200e+02,  -2.32881440e+02,   2.07022228e+03,
          6.07935200e+01,   1.44766400e+02,  -9.70688000e+00,
          8.69226000e+01,  -1.51670000e+02,   4.54800800e+01,
          4.59763024e+03,   9.35103080e+02],
       [  1.84748704e+03,   1.28695656e+03,  -8.56694920e+02,
         -1.27599448e+03,   1.62282716e+03,  -7.54124160e+02,
         -7.51028000e+00,   4.01664000e+01,  -3.29699200e+01,
          1.14110232e+03,  -1.35687120e+02,  -2.63592000e+00,
          6.56009360e+02,   8.72154800e+02,   1.00215168e+03,
          8.66088000e+00,  -2.59326412e+03,  -1.01746512e+03,
          7.01029200e+01,   1.19599640e+02,  -1.11451300e+03,
         -3.28862400e+01,  -7.56048800e+01,   2.46856000e+00,
         -4.24676000e+01,   8.73200800e+01,  -2.64428800e+01,
         -2.52171772e+03,  -4.01915040e+02],
       [  6.05564964e+03,   3.98362824e+03,  -2.65181920e+03,
         -3.94959140e+03,   4.86908816e+03,  -2.31025836e+03,
         -2.50621600e+01,   1.20792080e+02,  -1.27423720e+02,
          3.27100936e+03,  -4.23860120e+02,   3.68192000e+00,
          2.03039060e+03,   2.62953940e+03,   3.04463404e+03,
          7.71320400e+01,  -7.82742720e+03,  -3.14948508e+03,
          2.17003160e+02,   3.88421640e+02,  -3.43761624e+03,
         -1.04746440e+02,  -2.33571800e+02,   3.34929200e+01,
         -1.54264080e+02,   2.60516760e+02,  -8.35544800e+01,
         -7.80518924e+03,  -1.46023692e+03],
       [ -2.98017952e+03,  -2.07183312e+03,   1.37917192e+03,
          2.05426032e+03,  -2.60939344e+03,   1.21352736e+03,
          1.21126800e+01,  -6.46009600e+01,   5.35133600e+01,
         -1.83255016e+03,   2.18509400e+02,   4.05848000e+00,
         -1.05622988e+03,  -1.40291612e+03,  -1.61230440e+03,
         -1.48113600e+01,   4.17107144e+03,   1.63795232e+03,
         -1.12842480e+02,  -1.92861480e+02,   1.79391092e+03,
          5.29694400e+01,   1.21712560e+02,  -4.39320000e+00,
          6.87431200e+01,  -1.40352280e+02,   4.25722000e+01,
          4.06002808e+03,   6.50863040e+02]])
        
        _eps = 5.0e-8
        print 'Energy deviation:', numpy.linalg.norm(Energy - savedEnergy)
        print 'Force deviation:', numpy.linalg.norm(Force - savedForce)
        print 'Energy derivative deviation:', numpy.linalg.norm(EnergyDeriv -\
            savedEnergyDeriv)
        print 'Force derivative deviation:', numpy.linalg.norm(ForceDeriv -\
            savedForceDeriv)
        print 'Test tolerance:', _eps
        
        self.assertTrue(abs(Energy - savedEnergy) < _eps)
        self.assertTrue(numpy.linalg.norm(Force - savedForce) < _eps)
        self.assertTrue(numpy.linalg.norm(EnergyDeriv - savedEnergyDeriv) \
            < _eps)
        self.assertTrue(numpy.linalg.norm(ForceDeriv - savedForceDeriv) < _eps)

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        Output = MOPACOutput(sys.argv[1])
        Output.CalculateSpinPopulations()
    else:
        unittest.main()

