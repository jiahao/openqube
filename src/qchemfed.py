#!/usr/bin/env python

"""
Parses Q-Chem output files with FED calculations

For FED calculations, we grep TDDFT calculations
for the energies and transition dipoles of bright
excited states
"""

from numpy import array, dot, zeros
from numpy.linalg import norm, svd

StrengthThreshold = 0.1
eV = 0.03674932534 #converts eV to Hartree
Angstrom = 1.889725989
Debye = 1/2.5417462

SpinMultiplicity = {
    "Singlet":0,
    "Triplet":2,
}

class ElectronicState:
    def __init__(self):
        self.Amplitudes = {} # Ideally we want this to be a sparse matrix!
                             # But we will make do with a dictionary[(OCC, VIRT)] = ampl
        self.ExcitationEnergy = None
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

class ModelChemistry:
    def __init__(self):
        self.Exchange = None
        self.Correlation = None
        self.Basis = None
        self.Excited = None
        
    def __repr__(self):
        method = '' if self.Excited is None else self.Excited+'-'
        method += '' if self.Exchange is None else self.Exchange
        method += '' if self.Correlation is None else self.Correlation
        basis = '?' if self.Basis is None else self.Basis
        if method+basis=='?':
            return ''
        else:
            return method+'/'+basis

def MomentOfInertiaTensor(R, Weights = None, Center = True):
    "Moment of inertia tensor for a bunch of points, optionally weighted"
    I = zeros((3,3))
    if Center == True:
        C = Centroid(R, Weights)
    elif Center is None:
        C = zeros(3)
    else:
        C = Center

    for idx, p in enumerate(R):
        weight = 1 if Weights is None else Weights[idx]
        for i in range(3):
            for j in range(3):
                I[i,j] += weight * (p[i]-C[i]) * (p[j]-C[j])
    return I

def Centroid(R, Weights = None):
    "Centroid for a bunch of points, optionally weighted"
    C = zeros((3,))
    
    for idx, p in enumerate(R):
        weight = 1 if Weights is None else Weights[idx]
        C += weight * p
 
    normalization = len(R) if Weights is None else sum(Weights)
    if abs(normalization) < 0.001:
        normalization = 1.0

    return C/normalization

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
        for line in open(filename):
            if skiplines > 0:
                skiplines -= 1
                continue

            if mode == 'seek':
                if 'Standard Nuclear Orientation (Angstroms)' in line:
                    mode, skiplines = 'geometry', 2
                elif 'TDDFT Excitation Energies' in line:
                    mode, skiplines = 'tddft', 2
                elif 'Charges' in line:
                    mode, skiplines = 'charges', 3
                    Charges = []

                #Read in model chemistry
                elif len(line.split()) > 1 and line.split()[0].lower() == 'rpa' and line.split()[1].lower() == 'true':
                    QCModel.Excited = 'TD'
                elif len(line.split()) > 1 and line.split()[0].lower() == 'exchange':
                    QCModel.Exchange = line.split()[1].upper()
                elif len(line.split()) > 1 and line.split()[0].lower() == 'correlation':
                    QCModel.Correlation = line.split()[1].upper()
                elif len(line.split()) > 1 and line.split()[0].lower() == 'basis':
                    QCModel.Basis = line.split()[1].upper()

                #To trigger FED requires a two-stage trigger
                elif 'In RPA/TDDFT Excited States:' in line:
                #elif 'In TDDFT Excited States:' in line:
                    fed_trigger = True
                elif fed_trigger is True and 'FED Couplings Between Singlet Excited States' in line:
                    fed_trigger = False
                    mode = 'fed'
                    skiplines = 3

                    #Initialize
                    NumStates = len(States)
                    FEDCoupling = zeros((NumStates, NumStates))

            elif mode == 'geometry':
                if '----------------------------------------------------' in line:
                    mode = 'seek'
                    continue

                t = line.split()
                Geometry.append([t[1]] + list(map(lambda x:float(x)*Angstrom, t[2:5])))
            elif mode == 'tddft':
                t = line.split()
                if 'Excited state' in line and 'excitation energy (eV)' in line:
                    # Parse line of the form
                    # Excited state   2: excitation energy (eV) =    2.0813
                    State = ElectronicState()
                    State.Index, State.ExcitationEnergy = int(t[2][:-1]), float(t[-1])*eV
                elif 'Total energy for state' in line:
                    # Parse line of the form
                    #    Total energy for state   2:  -3332.812081058783
                    State.Energy = float(t[-1])
                elif 'Trans. Mom.:' in line:
                    # Parse line of the form
                    #    Trans. Mom.:  0.0371 X   0.0059 Y  -0.0349 Z
                    State.TransitionDipole = array(map(float, (t[2], t[4], t[6]))) 
                elif len(t)>0 and t[0] == 'Multiplicity':
                    # Parse line of the form
                    #    Multiplicity: Singlet
                    State.Multiplicity = SpinMultiplicity[t[-1]]
                elif len(t)>0 and t[0] == 'Strength':
                    # Parse line of the form
                    #     Strength   :  0.8578
                    State.Strength = float(t[-1])
                elif len(t)>2 and t[0] == 'X:' and t[2] == '-->':
                    # Parse line of the form
                    #    X: D(265) --> V(  1) amplitude = -0.2404
                    amplitude = float(t[-1])
                    mo_occ = int(t[1][2:5])
                    try: #In the event that we get a 3 digit orbital index
                        mo_virt = int(t[3][2:5])
                    except ValueError: #Will be parsed as separate token
                        mo_virt = int(t[4][:-1])

                    State.Amplitudes[mo_occ, mo_virt] = amplitude
                elif len(t) == 0:
                    # Terminate
                    States.append(State)
                elif '---------------------------------------------------' in line:
                    mode = 'seek'

            elif mode == 'fed':
                # The header for this section is
                #------------------------------------------------------------------------------
                #  States        X12(D)      X12(A)        dX12    Coupling(eV)
                #------------------------------------------------------------------------------
                #
                # Lines have the form:
                #
                #  1    2     -0.000001    0.000002   -0.000003   -0.3642343E-06
              
                # First check for termination
                if '------------------------------------------------------------------------------' in line:
                    mode = 'seek'
                    continue

                t = line.split()
                state1, state2 = int(t[0]) - 1, int(t[1]) -1
                FEDCoupling[state1, state2] = float(t[-1])*eV

            elif mode == 'charges':
                if '----------------------------------------' in line:
                    mode = 'seek'
                    continue
                
                t = line.split()
                Charges.append(float(t[-1]))

            else:
                raise ValueError, 'Encountered unknown mode '+mode
        #Pickle
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



def nearest(x, things):
    y = sorted([(abs(z-x), z) for z in things])
    return y[0][1]



def ForsterCoupling(d1, d2, r):
    """
    This function calculates the Forster coupling between two chromophores
    using the dipole-dipole coupling approximation.

    @param d1 Transition dipole of 1st chromophore in atomic units (3-vector)
    @param d2 Transition dipole of 2nd chromophore in atomic units (3-vector)
    @param r  Displacement vector between the two chromophores in atomic units
    @returns The coupling matrix element in atomic units
    @f[
    V = \frac {3 (d_1 \cdot \hat r) (d_2 \cdot \hat r) - (d_1 \cdot d_2)} {\vert r \vert^3 }
    @f]

    @note the formula doesn't care which direction the displacemnent vector is in
    """
    normr = norm(r)
    rn = r / normr ##Normalized distance vector
    Coupling = (3 * dot(d1, rn) * dot(d2, rn) - dot(d1, d2)) / normr**3
    return Coupling

def ForsterOrientationFactor(d1, d2, r):
    """
    This function calculates the Forster orientation factor between two chromophores
    using the dipole-dipole coupling approximation.

    @param d1 Transition dipole of 1st chromophore in atomic units (3-vector)
    @param d2 Transition dipole of 2nd chromophore in atomic units (3-vector)
    @param r  Displacement vector between the two chromophores in atomic units
    @returns The coupling matrix element in atomic units
    @f[
    \kappa
    @f]

    @note the formula doesn't care which direction the displacemnent vector is in
    """
    rn =  r / norm(r) ##Normalized distance vector
    d1n = d1/ norm(d1)
    d2n = d2/ norm(d2)
    Factor = 3 * dot(d1n, rn) * dot(d2n, rn) - dot(d1n, d2n)
    return Factor

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


