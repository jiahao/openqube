#!/usr/bin/env python
"""
parseci.py: Groks CI vector and orbitals from MOLPRO output file
Jiahao Chen, 20061011
"""

import sys
eV = 0.0367493245
Debye = 0.3934302

#---Parameters controlling program output---#000000#FFFFFF---------------------
THRESH_AO = 0.05 ** 0.5 #Threshold for ignoring atomic orbital coefficients
THRESH_CI = 0.05 ** 0.5 #Threshold for ignoring CI vector coefficients

THRESH_DM = 0.01*Debye #Threshold for ignoring transition dipole moment (in au)

DEBUG = False
ECHOINPUT = False
VERBOSE = True

#Converts MOLPRO configuration string into occupation number representation dictionary
def onr(molpro_config):
    occnos = {'alpha':[],'beta':[]}
    for orbital in molpro_config:
        if orbital=='2':
            occnos['alpha'].append(1)
            occnos['beta'].append(1)
        elif orbital in ['a', '+', '/']:
            occnos['alpha'].append(1)
            occnos['beta'].append(0)
        elif orbital in ['b', '-', '\\']:
            occnos['alpha'].append(0)
            occnos['beta'].append(1)
        elif orbital=='0':
            occnos['alpha'].append(0)
            occnos['beta'].append(0)
        else: #String is badly formatted, raise an exception!
            raise ValueError
    return occnos

#Calculates magnetic spin number of a given configuration
def spin_number(molpro_config):
    totalspin = 0
    for orbitalocc in molpro_config:
        if orbitalocc == 'a':
            totalspin += 0.5
        elif orbitalocc == 'b':
            totalspin -= 0.5
    return totalspin

#Counts number of electrons in a given configuration
def numelecs(molpro_config):
    num = 0
    for orbitalocc in molpro_config:
        if orbitalocc == 'a':
            num += 1
        elif orbitalocc == 'b':
            num += 1
        elif orbitalocc == '2':
            num += 2
    return num

#Converts a difference between two ONR configuration strings into orbital (de)population notation
def deltaorbs(refconfig, baseconfig):
    if len(baseconfig) != len(refconfig): #Must be same length!
        raise IndexError
    
    ref_onr = onr(refconfig)
    exc_onr = onr(baseconfig)

    #determine nature of excitations
    holes = []
    particles = []
    for ii in range(len(refconfig)):
        for spin in ['alpha', 'beta']:
            popchange = exc_onr[spin][ii]-ref_onr[spin][ii]
            if popchange == 1: #Found a particle-type excitation
                particles.append((ii,spin))
            elif popchange == -1: #Found a hole-type excitation
                holes.append((ii,spin))

    output = '('
    for hole in holes:
        output+='%d%s ' % (hole[0]+1, hole[1][0])
    output = output[:-1] + ') --> ('
    for particle in particles:
        output+='%d%s ' % (particle[0]+1, particle[1][0])
    output = output[:-1] + ')'
    return output

#---Parse command-line options---
for options in sys.argv: #command line
    if options in ['-d', '--debug']:
        DEBUG = True
        sys.argv.pop(0)
    elif options in ['-e', '--echo']:
        ECHOINPUT = True
        sys.argv.pop(0)
    elif options in ['-v', '--verbose']:
        VERBOSE = True
        sys.argv.pop(0)
    elif options in ['-q', '--quiet']:
        VERBOSE = False
        ECHONINPUT = False
        DEBUG = False
        sys.argv.pop(0)
    elif options in ['-h,' '--help', '-?']:
        print """parseci.py : Simplifies MOLPRO output
Options
-d, --debug   : Print debugging information
-e, --echo    : Echo MOLPRO output file as it is being read
-v, --verbose : Print extra information (default)
-q, --quiet   : Minimal information
-h, --help    : Print this help screen
    """
    exit()

#---Get file name and open file---#000000#FFFFFF-------------------------------
try:
    mpoutfn=sys.argv[1]
except IndexError:
    mpoutfn='out.out'

print 'Opening %s' % mpoutfn
mpoutf=open(mpoutfn,'r')

for line in mpoutf:
    if ECHOINPUT:
        print 'Read: ',line,
    toks = line.split()
    #---Read XYZ geometry---#000000#FFFFFF-------------------------------------
    if 'ATOMIC COORDINATES' in line:
        if VERBOSE:
            print 'Reading geometry'
        #Skip to start of XYZ geometry
        for _ in range(4):
            line = mpoutf.next()
            if ECHOINPUT:
                print 'Read: ',line,

        #---Read number of atoms---
        #numatoms = int(line)
        #if VERBOSE:
        #    print 'Geometry has %d atoms' % numatoms
        
        #---Read coordinates---
        atoms = []
        xyz = []
        numatoms = 0
        #If there is a comment line, skip it
        #comment  = mpoutf.next().strip()
        #toks = comment.replace(',',' ').split()
        toks = line.split()
        try: #Test if there are three numbers following an atom
            coord = map(float, toks[3:5])
        except (IndexError, ValueError): #Oops, doesn't work, it's a comment
            if ECHOINPUT:
                print 'Read: ',line,
            if VERBOSE:
                print line
            start=0        
        else: #Format matches xyz
            while len(toks)>0: #Read until line is blank
                numatoms+=1
                if DEBUG:
                    print toks
                atoms.append(toks[1])
                xyz.append(coord)
                line=mpoutf.next()
                if ECHOINPUT:
                    print 'Read: ',line,
                toks=line.split()
        if VERBOSE:
            print 'Finished reading geometry with %d atoms' % numatoms
    #---New program run, initialize state data---#000000#FFFFFF----------------
    if 'PROGRAM' == line[1:8]:
        numstates=0
        state={}
        statenrgs={}
        print 80*'-'
        print 'New %s program output section' % toks[2]
    
    #---Read energies---#000000#FFFFFF-----------------------------------------
    if len(toks)>3 and '!'==line[1] and 'ENERGY'==toks[-2]:
        numstates += 1
        statenum = int(toks[-3].replace('STATE','')[:-2])
        statenrg = float(toks[-1])
        statenrgs[statenum]=statenrg
        print 'State %2d has energy %f' % (statenum, statenrg)

    #---Read transition dipole moments---#000000#FFFFFF------------------------
    if len(toks)>3 and toks[1]=='trans' in line:
        #Python chokes on < and >, redo tokenize
        toks=line.replace('<','').replace('>','').split()

        #---Break up printed matrix element---
        matelt = toks[2].split('|')
        fromstate = int(matelt[2].split('.')[0])
        tostate = int(matelt[0].split('.')[0])
        direction = matelt[1][-1:]
        dipole = float(toks[3])
        if abs(dipole) > THRESH_DM:
            print 'Transition dipole from state %2d to state %2d in %s direction is %f D' % (fromstate, tostate, direction, dipole/Debye) 
        
    #---Read in orbitals---#000000#FFFFFF-------------------------------------
    if 'ORBITALS' in toks:
        #These variables track the estimated beginning and end of active space
        orbitalactivespacebegin=0
        orbitalactivespaceend=0
        #---Read in type of orbitals---
        orbtype = toks[0] #format for ELECTRON ORBITALS and NATURAL ORBITALS differ
        if VERBOSE:
            print 'Found %s orbitals' % orbtype
        #Skip ahead to header
        while not 'Orb' in line:
            line = mpoutf.next()
            if ECHOINPUT:
                print 'Read: ',line,
        mpoutf.next() #Skip one more line
        if ECHOINPUT:
            print 'Read: ',line,
        line=mpoutf.next()
        if ECHOINPUT:
            print 'Read: ',line,
        toks=line.split()
        #---Read in AO key---
        aokey=[]
        while len(toks)>0:
            while len(toks)>0:
                aokey.append((int(toks.pop(0)), toks.pop(0)))
                
            line=mpoutf.next()
            if ECHOINPUT:
                print 'Read: ',line,
            toks=line.split()
        if DEBUG:
            print 'aokey=',aokey

        #Skip ahead to MO specs
        while len(toks)==0:
            line = mpoutf.next()
            if ECHOINPUT:
                print 'Read: ',line,
            toks = line.split()

        #---Read AO coefficients---
        mo = {}
       
        while len(toks)>0: #Loop over molecular orbital
            try:
                orbnum=int(toks.pop(0).split('.')[0])
            except ValueError:
                break
            #Read occupation number
            tok = toks.pop(0)
            #In old MOLPRO versions, sometimes have '+' and '-' occupations
            if tok == '-':
                occno = 1
            elif tok == '+':
                occno = 1
            else:
                occno =float(tok)
            if occno == 2: #Population is not exactly 2, guess start of active space
                orbitalactivespacebegin = orbnum
            if occno > 0: #Population is not exactly 0, guess end of active space
                orbitalactivespaceend = orbnum
            orbnrg=float(toks.pop(0).replace('\x00',''))
            if orbtype == 'ELECTRON': #Has entry for Coulson energy
                orbcoul=float(toks.pop(0))

            aoid = 0
            mo[orbnum]=[]
            thismo=mo[orbnum]
            if VERBOSE:
                print 'Orbital %d has occupation %f and energy %f' % (orbnum, occno, orbnrg)
            while len(toks)>0: #Loop over line entry
                while len(toks)>0: #Loop over entry in line
                    aoid+=1
                    #Sometimes split() gives bizarre null bytes
                    tok = toks.pop(0).replace('\x00','')
                    while tok == '':
                        tok = toks.pop(0).replace('\x00','')
                    coeff=float(tok)
                    if abs(coeff) > THRESH_AO:
                        thismo.append((aoid,coeff))
                        print 'Orbital %d has AO (#%2d: %s on atom %d) with coefficient %f (weight %f)' % (orbnum, aoid, aokey[aoid][1], aokey[aoid][0], coeff, coeff**2)
                line=mpoutf.next()
                if ECHOINPUT:
                    print 'Read: ',line,
                toks=line.split()
            line=mpoutf.next()
            if ECHOINPUT:
                print 'Read: ',line,
            toks=line.split()

        #---Print estimated active space---
        orbitalactivespacebegin+=1
        print 'Lowest orbital with occupation less than 2 = %d' % orbitalactivespacebegin
        print 'Highest orbital with occupation greater than 0 = %d' % orbitalactivespaceend
        
        if VERBOSE:
            print 'Finished reading orbitals'

    #---Read in states---#000000#FFFFFF----------------------------------------
    if 'CI vector' in line:
        if VERBOSE:
            print 'Found CI vector'
        #Initialize state data
        for i in range(0,numstates):
            state[i]=[] #state is a list of lists
        #Skip two lines, process the thirs
        for _ in range(3):
            line=mpoutf.next()
            if ECHOINPUT:
                print 'Read: ',line,
        toks=line.split()
        #---Read CI vector---
        while len(toks)>0: #Iterate over configurations
            config=toks.pop(0)
            for i in range(0,numstates):
                #Advance to start of configs and coefficients
                if len(toks)==0:
                    line=mpoutf.next()
                    if ECHOINPUT:
                        print 'Read: ',line,
                    toks=line.split()

                try: #hack to not crash when non-unique state numbers exist
                    coeff=float(toks.pop(0))
                except ValueError:
                    print 'WARNING: Found multiple states with same state numbers!!'
                    numstates = i
                else:    
                    if abs(coeff) > THRESH_CI:
                        state[i].append((config,coeff))
                        if VERBOSE:
                            print "State %2d has weight %f for configuration %s" % (i+1, coeff, config)

            line=mpoutf.next()
            if ECHOINPUT:
                print 'Read: ',line,
            toks=line.split()

        #---Read state energies---
        line=mpoutf.next()
        if ECHOINPUT:
            print 'Read: ',line,
        toks=line.split()
        #Skip text label
        toks.pop(0)
        toks.pop(0)
        for i in range(0,numstates):
            if len(toks)==0:
                line=mpoutf.next()
                if ECHOINPUT:
                    print 'Read: ',line,
                toks=line.split()
            statenrg=float(toks.pop(0))
            statenrgs[i]=statenrg
            if VERBOSE:
                print "State %2d has energy %f" % (i+1, statenrg)

        if VERBOSE:
            print 'Finished reading CI vector'

        #---Identify dominant ground state configuration---
        #Assumed only ONE dominant configuration
        #This may have to be changed!
        maxgsconfig = ''
        maxgsweight = 0
        for config in state[0]:
            if abs(config[1]) > abs(maxgsweight):
                maxgsconfig = config[0]
                maxgsweight = config[1]

        print 'The ground state has dominant configuration %s with coefficient %f (probability %f)' % (maxgsconfig, maxgsweight, maxgsweight**2)
        if abs(maxgsweight) < 0.5**0.5:
            print 'Warning, ground state appears to not have any one dominant configuration!'
            if VERBOSE:
                print 'The configurations with weights > %f are' % THRESH_CI
                for config in state[0]:
                    print '\tCoefficient %f (Density %f) \tConfiguration %s' % (config[1], config[1]**2, config[0])

        #Convert configuration into more convenient data structure in occupation number representation
        
        #---Print results of state analyses---
        numorbs = orbitalactivespaceend - orbitalactivespacebegin - 1
        print 'The effective active space is %d electrons in %d orbitals' % (numelecs(maxgsconfig), numorbs)
        print 'Excited states:'
        for i in range(1,numstates):
            print 'State %d:\tExcitation energy: %f eV' % (i+1,(statenrgs[i] - statenrgs[0])/eV)
            #Convert ONR configuration into orbital (de)population notation:
            for config in state[i]:
                excitations=deltaorbs(maxgsconfig, config[0])
                print '\tCoefficient %f (Density %f) \tExcitation %s\tm_s %d --> %d' % (config[1], config[1]**2, excitations, spin_number(maxgsconfig), spin_number(config[0]))
            #Compute norm
            norm = 0
            for config in state[i]:
                norm += config[1]**2
            print 'Norm: %f' % norm
            print
    #---Detect any errors in output---#000000#FFFFFF---------------------------
    if 'ERROR EXIT' in line[1:11]:
        print 'Error detected in MOLPRO output file!'
        
mpoutf.close()
