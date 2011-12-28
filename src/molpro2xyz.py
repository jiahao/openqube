#!/usr/bin/env python
print """
optg2xyz.py: Extracts result of MOLPRO's OPTG program and outputs a XYZ trajectory
Jiahao Chen, 20060906
"""
import sys

try:
    outname=sys.argv[1]
except IndexError:
    print 'Error: did not receive a output filename'
else:
    try:
        xyzname=sys.argv[2]
    except IndexError:
        xyzname=outname+'-optg.xyz' #Generate default output file name
        print 'No output file name specified, using default'
print 'Writing output to file '+xyzname

try: #load MolPro output file
    outfile=open(outname+'.out','r')
except IOError:
    print 'FATAL ERROR: Could not load output file %s.out' % outname
else:
    print 'Processing %s.out' % outname
    skip=0
    coords=[]
    atoms=[]
    numatoms=''
    comment=''
    xyztraj=''
    readmode=0
    numgeoms=1
    for line in outfile:
        if skip == 0:
            toks = line.replace(',',' ').split()

            #Search for initial geometry
            if 'geomtyp=xyz' in toks:
                print 'Initial geometry specified in XYZ format'
                readmode=10
            elif 'geometry={' in toks:
                if readmode==10:
                    readmode=1
                else:
                    print "I don't understand the input, not sure if initial geometry is in XYZ format"
                    print "Skipping initial geometry"

            #Reading XYZ geometry
            if readmode==1: #Read number of atoms
                try:
                    numatoms=int(line)
                except ValueError:
                    "" #Do nothing, this is the wrong line
                else:
                    xyztraj = xyztraj + str(numatoms)
                    readmode=2
            elif readmode==2: #Read comment line
                comment=line[:-1]
                xyztraj = xyztraj + '\n' + comment
                readmode=3
            elif readmode==3: #Read coordinates
                if '}' in toks:
                    readmode=0 #End of geometry read mode
                else:
                    atoms.append(toks[0])
                    xpos = toks[1]
                    ypos = toks[2]
                    zpos = toks[3]

                    xyztraj = xyztraj + '\n%s\t%s\t%s\t%s' % (toks[0],xpos, ypos, zpos)

                    
            #Reading OPTG
            if 'cartesian' in toks: #Used cartesian coordinates
                print 'OPTG used cartesian coordinates'
                readmode=11
            if 'Optimization' in toks: # Start of a geometry
                if readmode==11:
                    skip=3
                else:
                    print 'OPTG did not specify cartesian coordinates'
                    print 'I only understand cartesian coordinates!'
                    skip=9999999999
            if 'ANGSTROM' in toks: #Read coordinate
                coords.append(toks[5])
            if 'Convergence:' in toks: #finished reading out a geometry
                numgeoms=numgeoms+1
                if coords!=[]: #Have geometry, write coordinates
                    xyztraj=xyztraj+'\n%d\n%s' % (numatoms, comment)
                    for atom in atoms:
                        xpos = coords.pop(0)
                        ypos = coords.pop(0)
                        zpos = coords.pop(0)
                        xyztraj = xyztraj + '\n%s\t%s\t%s\t%s' % (atom, xpos, ypos, zpos)
            if 'ERROR' in toks:
                print 'Warning, MOLPRO error detected in output'
        else: #Skipped line processing
            skip = skip - 1
    outfile.close()

    #print xyztraj
    print 'Writing trajectory of %d geometries (including original) to %s' % (numgeoms, xyzname)
    xyzfile=open(xyzname,'w')
    xyzfile.write(xyztraj+'\n')
    xyzfile.close
