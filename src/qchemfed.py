#!/usr/bin/env python

"""
Parses Q-Chem output files with FED calculations

For FED calculations, we grep TDDFT calculations
for the energies and transition dipoles of bright
excited states
"""

from numpy import array, dot, zeros
from numpy.linalg import norm, svd
from data import eV, Angstrom
from analysis_geometry import Centroid, MomentOfInertiaTensor
from analysis_exciton import SpinMultiplicity, nearest, ForsterCoupling, ForsterOrientationFactor

StrengthThreshold = 0.1

def ParseFED(filename, doPlot = False):
    """
    See Q-Chem 3.2 manual,
    10.17 Electronic Couplings for Electron Transfer and Energy Transfer
    """

    #Pickle
    try:
        import pickle
        f = open(filename+'.pkl', 'rb')
        Geometry = pickle.load(f)
        States = pickle.load(f)
        QCModel = pickle.load(f)
        Charges = pickle.load(f)
        f.close()
        print 'Loaded from pickle', filename
    except (EOFError, IOError):
        mode = 'seek'
        skiplines = 0
        fed_trigger = False
        #Molecular geometry
        Geometry = []
        #List of states
        States = []
        #FED coupling matrix; will be initialized later
        FEDCoupling = None
        QCModel = ModelChemistry()
        
        f = open(filename+'.pkl', 'wb')
        pickle.dump(Geometry, f)
        pickle.dump(States, f)
        pickle.dump(QCModel, f)
        f.close()
        #print 'Cached parsed output in', filename+'.pkl'

    ###################
    # Post-processing #
    ###################
    if doPlot:
        print 'Perceived the following state information:'
        for state in States:
            print state
            state.isValid()

    Coords = array([l[1:] for l in Geometry])

    if len(States) == 0:# or FEDCoupling is None:
        print filename, 'No data found'
        return False

    if FEDCoupling is None: #Assume this is a monomer calculation
        if doPlot:
            import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d import Axes3D

            print 'Plotting geometry in Figure 1'
            fig = plt.figure(1)
            ax = Axes3D(fig)
            for x1, y1, z1 in Coords/Angstrom:
                for x2, y2, z2 in Coords/Angstrom:
                    if 0 < ((x1-x2)**2 +  (y1-y2)**2 + (z1-z2)**2)**0.5 < 1.5:
                        ax.plot((x1,x2), (y1,y2), (z1,z2), 'k-', marker='.')
            plt.title('Geometry in Angstroms')

            c = Centroid(Coords)/Angstrom
            for state in States:
                tdip = state.TransitionDipole/Angstrom

                if norm(tdip) > 0.2:
                    d = c + tdip
                    ax.plot((c[0], d[0]), (c[1], d[1]), (c[2], d[2]))
                    ax.text(d[0], d[1], d[2], str(state.Index))

            print 'Plotting absorption spectrum in Figure 2'
            fig = plt.figure(2)
            HartreeToNm = 1239.8/27.2113839
            x = [HartreeToNm/state.ExcitationEnergy for state in States]
            y = [state.Strength for state in States]
            plt.stem(x, y, markerfmt='*', linefmt='k-')
            for idx, state in enumerate(States):
                plt.text( x[idx], y[idx], str(state.Index), horizontalalignment='left')
                
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('Oscillator strength')
            plt.title('Absorption spectrum predicted from '+str(QCModel))

            plt.show()

        Inertia = MomentOfInertiaTensor(Coords)

        #Save monomer data
        f = open('monomer.pkl', 'wb')
        pickle.dump(Inertia, f)
        pickle.dump(len(States), f)
        for state in States:
            pickle.dump(state.TransitionDipole, f)
        f.close()

        print 'Saved monomer data in monomer.pkl:'
        print 'Moment of inertia tensor:\n', Inertia
        print 'Transition dipoles:'
        for state in States:
            print state.Index, state.TransitionDipole
        return False
    else:
        #Load data from monomer pickle
        try:
            f = open('monomer.pkl', 'rb')
            #print 'Loading monomer data from monomer.pkl:'
            MonomerInertia = pickle.load(f)
            #print 'Moment of inertia tensor:\n', MonomerInertia
            numstates = pickle.load(f)
            MonomerTransitionDipoles = []
            #print 'Transition dipoles:'
            for _ in range(numstates):
                tdip = pickle.load(f)
                MonomerTransitionDipoles.append(tdip)
                #print i, tdip
            f.close()

        except (EOFError, IOError):
            raise ValueError, 'Tried to analyze FED calculation without monomer data in monomer.pkl.'

    #Assume DIMER calculation
    # Step 1. Calculate geometric parameters

    # Find distance between centroids
    numatoms = len(Geometry)/2
    r1 = Centroid(Coords[:numatoms,:])
    r2 = Centroid(Coords[numatoms:,:])
    r = r1 - r2
    distance = norm(r)

    # Calculate local inertia tensors
    I1 = MomentOfInertiaTensor(Coords[:numatoms,:])
    I2 = MomentOfInertiaTensor(Coords[numatoms:,:])
    
    _, _, Vm = svd(MonomerInertia)
    U1, _, _ = svd(I1)
    U2, _, _ = svd(I2)
    M1 = dot(U1, Vm)
    M2 = dot(U2, Vm)
    if doPlot:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        print 'Plotting geometry in Figure 1'
        fig = plt.figure(1)
        ax = Axes3D(fig)
        for x1, y1, z1 in Coords/Angstrom:
            for x2, y2, z2 in Coords/Angstrom:
                if 0 < ((x1-x2)**2 +  (y1-y2)**2 + (z1-z2)**2)**0.5 < 1.5:
                    ax.plot((x1,x2), (y1,y2), (z1,z2), 'k-', marker='.')
        plt.title('Geometry in Angstroms')

        #Transition dipoles on first monomer
        c = r1/Angstrom
        for idx, td in enumerate(MonomerTransitionDipoles):
            tdip = dot(M1, td/Angstrom)
            if norm(tdip) > 0.2:
                d = c + tdip
                ax.plot((c[0], d[0]), (c[1], d[1]), (c[2], d[2]), 'b-')
                ax.text(d[0], d[1], d[2], idx+1)
        
        #Transition dipoles on second monomer
        c = r2/Angstrom
        for idx, td in enumerate(MonomerTransitionDipoles):
            tdip = dot(M2, td/Angstrom)
            if norm(tdip) > 0.2:
                d = c + tdip
                ax.plot((c[0], d[0]), (c[1], d[1]), (c[2], d[2]), 'r-')
                ax.text(d[0], d[1], d[2], idx+1)

        #Transition dipoles on dimer
        c = Centroid(Coords)/Angstrom
        for state in States:
            tdip = state.TransitionDipole/Angstrom
            if norm(tdip) > 0.2:
                d = c + tdip
                ax.plot((c[0], d[0]), (c[1], d[1]), (c[2], d[2]), 'g-')
                ax.text(d[0], d[1], d[2], state.Index)

    #Calculate all possible Forster couplings
    ForsterCouplings = {}
    for idx, td in enumerate(MonomerTransitionDipoles):
        idx2, td2 = idx, td
        c = abs(ForsterCoupling(dot(M1, td), dot(M2, td2), r)/eV)
        ForsterCouplings[c] = (idx+1, idx2+1)
        #for idx2, td2 in enumerate(MonomerTransitionDipoles):
        #    if idx <= idx2:
        #        c = abs(ForsterCoupling(td, dot(M, td2), r)/eV)
        #        ForsterCouplings[c] = (idx+1, idx2+1)

    #Prune high energy states
    #States = filter(lambda x: x.ExcitationEnergy < 0.073, States)
    
    #Prune dark states
    States = filter(lambda x: norm(x.TransitionDipole) > 0.2, States)
    
    #Prune insignificant transition amplitudes
    for state in States:
        threshold = (len(state.Amplitudes)+1)**-0.5
        newamplitudes = dict()
        for (occ, virt), amplitude in state.Amplitudes.items():
            if abs(amplitude) > threshold:
                newamplitudes[(occ, virt)] = amplitude
        state.Amplitudes = newamplitudes

    #Match states
    MatchingStates = []
    for idx1, state1 in enumerate(States):
        for idx2, state2 in enumerate(States):
            if idx1 < idx2:
                MatchingStates.append((idx1, idx2))
            if len(MatchingStates) == 1:
                break
        if len(MatchingStates) == 1:
            break
   
    for idx1, idx2 in MatchingStates:
        TheFEDCoupling = FEDCoupling[idx1, idx2]
        print filename, 
        print distance/Angstrom, idx1+1, idx2+1,
        #print States[idx1].ExcitationEnergy/eV, States[idx2].ExcitationEnergy/eV,
        print abs(TheFEDCoupling)/eV,
        x =nearest(abs(TheFEDCoupling), ForsterCouplings.keys())
        print x/eV,
        n1, n2 = ForsterCouplings[x]
        print n1, n2, ForsterOrientationFactor(MonomerTransitionDipoles[n1], MonomerTransitionDipoles[n2], r,),
        print
    return States, FEDCoupling

if __name__ == '__main__':
    import glob, sys
    if len(sys.argv) == 1:
        print """I will now proceed to process *.out in this directory.
If you did not want this behavior, abort now and rerun
with specified output file(s) on command line"""
        args = ['*.out']
    else:
        args = sys.argv[1:]

    #Get glob of globs
    filenames = []
    for fileglob in args:
        filenames += glob.glob(fileglob)

    #Uniqfy
    filenames = dict().fromkeys(filenames).keys()
    filenames.sort()

    doThisPlot = True if len(filenames) == 1 else False
    for fname in filenames:
        data = ParseFED(fname, doThisPlot)


