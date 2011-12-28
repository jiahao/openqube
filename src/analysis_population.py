#!/usr/bin/env python
#######################################
# Functions for population analyses

import scipy
from MatrixIO import Inv, InvSqrt, Sqrt, PurifyReal, \
    MatrixFromSymmetricPartAsVector, SymmetricPartAsVector

from MatrixIO import PrintMatrix
from copy import deepcopy
from numpy.linalg.linalg import eigvals


_epsilon = scipy.sqrt(scipy.finfo(float).eps)

def StickyFragmentDensity(P, Px, S, Sx, Sxx, Fragment, Verbose = False):
    """Does the sticky electron calculation"""

    #InvSxx = Inv(Sxx)
    #From RI on AOx basis
    #M = 0.5 * Px * Sx * InvSxx
    #X = PI * P + M

    if Verbose: print 'And here begins the great quest for a solution'

    n = S.shape[0]

    ZZ = scipy.matrix(scipy.zeros((n * 2, n * 2)))
    ZZ[0:  n, 0:  n] = Px
    ZZ[n:2 * n, 0:  n] = PI * P
    ZZ[0:  n, n:2 * n] = P * PI
    ZZ[n:2 * n, n:2 * n] = 0

    SS = scipy.matrix(scipy.zeros((n * 2, n * 2)))
    SS[0:  n, 0:  n] = S
    SS[n:2 * n, 0:  n] = Sx
    SS[0:  n, n:2 * n] = -Sx
    SS[n:2 * n, n:2 * n] = Sxx

    #Try normalizing Sxx first
    Nm = scipy.matrix(scipy.zeros((n, n)))
    for i in range(n):
        Nm[i, i] = Sxx[i, i] ** -0.5

    NNm = scipy.matrix(scipy.eye(2 * n))
    NNm[n:, n:] = Nm

    # InvNNm = Inv(NNm)

    #ZZ = InvNNm * ZZ * InvNNm
    #SS = NNm * SS * NNm
    #end try

    SSS = deepcopy(SS)
    SSS[n:2 * n, 0:  n] = 0
    SSS[0:  n, n:2 * n] = 0

    SqrtSS = Sqrt(SSS)
    InvSqrtSS = Inv(SqrtSS)

    SSt = InvSqrtSS * SS * InvSqrtSS

    Sxt = SSt[0:n, n:]

    if Verbose:
        PrintMatrix(SS, 'SS')
        PrintMatrix(SSS, 'Transformation')
        PrintMatrix(SSt, 'Transformed SS (super-overlap)')

        PrintMatrix(Sxt * Sxt.T, 'Our approximation to the identity S0x * Sx0')
        PrintMatrix(Sxt.T * Sxt, 'Our approximation to the identity Sx0 * S0x')
        print 'Eigenvalues'
        print eigvals(Sxt * Sxt.T)
        PrintMatrix(Sx * Sx.T, 'The original product S0x * Sx0')
        print eigvals(Sx * Sx.T)

    ZZt = InvSqrtSS * ZZ * InvSqrtSS

    if Verbose:
        PrintMatrix(ZZ, 'ZZ')
        PrintMatrix(ZZt, 'Transformed ZZ (super density)')

    Pxt = deepcopy(ZZt[0:  n, 0:  n])

    ZZt[n:2 * n, 0:  n] += 0.5 * Pxt * Sxt
    ZZt[0:  n, n:2 * n] += 0.5 * Sxt.T * Pxt
    ZZt[0:n, 0:n] = 0

    if Verbose: PrintMatrix(ZZt, 'Split h&h')

    AA = SqrtSS * ZZt * SqrtSS

    if Verbose: PrintMatrix(AA, 'Fragment guess')

    X = AA[0:n, n:2 * n]

    return 0.5 * (X + X.T)



def StickyCharge(P, S, Px, Sx, Sxx, Fragment):
    """Does the sticky electron calculation"""
    Pf = StickyFragmentDensity(P, S, Px, Sx, Sxx, Fragment)
    return Fragment.Z - scipy.trace(Pf * S)



def TikhonovRegularizedMatrix(M, threshold = 0.0, Verbose = True):
    """Return a Tikhonov-regularized matrix."""

    if threshold > 0:
        regularize = threshold
    else:
        regularize = 0.0
        evals = eigvals(M)
        for eval in evals:
            if eval < 0:
                regularize = max(regularize, abs(eval))
        if regularize > 0.0:
            'Regularizing objective with threshold %f' % regularize
            regularize += _epsilon #Just to be safe

    return M + threshold * scipy.eye(M.shape[0])


#############
def StickyByFixedPointIteration(A, P, Px, S, Sx, Sxx, Fragment):
    PI = Fragment.ProjectionMatrix()

    #from numpy.linalg import pinv
    Sinv = Inv(S)

    B = S * Px * Sx + S * PI * P * Sxx - Sx * P * PI * Sx
    B = B + B.T

    B = Sinv * B * Sinv
    C = Sinv * Sxx

    from MatrixIO import Lyapunov

    #A = A * 0
    from Optimizer import DIISExtrapolator
    DIIS = DIISExtrapolator()
    DIIS.AddData(A*0, ObjectiveFunctionDerivative(A*0, P, Px, S, Sx, Sxx, Fragment))
    for i in range(1000):

        Ax = ObjectiveFunctionDerivative(A, P, Px, S, Sx, Sxx, Fragment)

        DIIS.AddData(A, Ax)

        Aex = DIIS.Extrapolate()
        Anorm = ObjectiveFunction(A, P, Px, S, Sx, Sxx, Fragment)
        Aexnorm = ObjectiveFunction(Aex, P, Px, S, Sx, Sxx, Fragment)

        if Aexnorm < Anorm: #Accept
            print Aexnorm, '<', Anorm, ': accept DIIS'
            A = scipy.reshape(Aex, A.shape)
            Anorm = Aexnorm
        else:
            print Aexnorm, '>', Anorm, ': reject DIIS'

        print 'Iteration', i, ':', Anorm
        
        D = 2 * Sinv * Sx * A * Sx * Sinv
        A = Lyapunov(C, B + D, CheckSane = False, Verbose = False)
    return A


##############################
# Functions that generate data

def ObjectiveFunction(A, P, Px, S, Sx, Sxx, Fragment):
    """Calculates objective function."""
    PI = Fragment.ProjectionMatrix()

    n = S.shape[0]
    SS = scipy.matrix(scipy.zeros((n * 2, n * 2)))

    SS[0:  n, 0:  n] = S
    SS[n:2 * n, 0:  n] = Sx
    SS[0:  n, n:2 * n] = -Sx
    SS[n:2 * n, n:2 * n] = Sxx

    ZZ = scipy.matrix(scipy.zeros((n * 2, n * 2)))

    #Here is a trick to reshape A in case it is passed in as a flat array
    AA = A.view()
    AA.shape = n, n

    ZZ[0:  n, 0:  n] = Px
    ZZ[n:2 * n, 0:  n] = PI * P - AA
    ZZ[0:  n, n:2 * n] = P * PI - AA
    ZZ[n:2 * n, n:2 * n] = 0

    SqrtSS = Sqrt(SS)
    YY = SqrtSS * ZZ * SqrtSS
    YYSq = YY.H * YY
    YYSq = TikhonovRegularizedMatrix(YYSq, threshold = _epsilon)
    try:
        assert scipy.trace(YYSq).imag < 1e-6
        return scipy.trace(YYSq).real
    except SyntaxError: #YYSq is real
        return scipy.trace(YYSq)

    #There are five possible solutions (not necessarily plausible!)
    #1
    #return scipy.trace(ZZ * SS * ZZ * SS)
    #2
    return scipy.trace(ZZ * SS * ZZ.T * SS)
    #3
    #return scipy.trace(ZZ * ZZ * SS * SS)
    #4
    #return scipy.trace(ZZ * ZZ.T * SS * SS)
    #5 - same as #4
    #return scipy.trace(ZZ.T * ZZ * SS * SS)



def ObjectiveFunctionDerivative(A, P, Px, S, Sx, Sxx, Fragment):
    """Calculates derivative of the objective function."""
    PI = Fragment.ProjectionMatrix()

    n = S.shape[0]
    SS = scipy.matrix(scipy.zeros((n * 2, n * 2)))

    SS[0:  n, 0:  n] = S
    SS[n:2 * n, 0:  n] = Sx
    SS[0:  n, n:2 * n] = -Sx
    SS[n:2 * n, n:2 * n] = Sxx


    ZZ = scipy.matrix(scipy.zeros((n * 2, n * 2)))

    #Here is a trick to reshape A in case it is passed in as a flat array
    AA = A.view()
    AA.shape = n, n
    ZZ[0:  n, 0:  n] = Px
    ZZ[n:2 * n, 0:  n] = PI * P - AA
    ZZ[0:  n, n:2 * n] = P * PI - AA
    ZZ[n:2 * n, n:2 * n] = 0

    SqrtSS = Sqrt(SS)
    YY = SqrtSS * ZZ * SqrtSS
    #YYSq = YY.H * YY
    X = -2 * SqrtSS * YY.conj() * SqrtSS

    #1 tr (ZZ * SS * ZZ * SS)
    #X = -2 * SS * ZZ.T * SS
    #2 tr (ZZ * SS * ZZ.T * SS)
    #X = -2 * (SS * ZZ * SS)
    #3 tr (ZZ * ZZ * SS * SS)
    #X = -(ZZ.T * SS * SS + SS * SS * ZZ.T)
    #4 tr (ZZ * ZZ.T * SS * SS)
    #X = -2 * SS * SS * ZZ
    #5 tr(ZZ.T * ZZ * SS * SS)
    #X = -2 * ZZ * SS * SS

    #PrintMatrix(X, 'X')
    Y = X[0:n, n:2 * n] + X[n:2 * n, 0:n]
    return Y



###############################################
# Helper functions for optimizers

def ObjectiveFunctionSymmetricOnly(A, P, Px, S, Sx, Sxx, Fragment):
    """Identical to ObjectiveFunction, but A is input differently
    Here A is a vector containing one half of the elements. We
    assume that the fragment density is symmetric"""
    AA = MatrixFromSymmetricPartAsVector(A)
    return ObjectiveFunction(AA, P, Px, S, Sx, Sxx, Fragment)



def ObjectiveFunctionDerivativeSymmetricOnly(A, P, Px, S, Sx, Sxx, Fragment):
    """Identical to ObjectiveFunctionDerivative, but A is input differently
    Here A is a vector containing one half of the elements. We
    assume that the fragment density is symmetric
    
    Also returns symmetric part of derivative only"""
    AA = MatrixFromSymmetricPartAsVector(A)
    Y = ObjectiveFunctionDerivative(AA, P, Px, S, Sx, Sxx, Fragment)
    #PrintMatrix(Y, 'Y')
    for i in range(AA.shape[0]):
        for j in range(AA.shape[1]):
            if i != j:
                Y[i, j] *= 2
    #PrintMatrix(Y, 'Y')
    return SymmetricPartAsVector(Y)



def LowdinCharge(P, S, Fragment):
    """Calculates Lowdin charge for the fragment"""
    M = Sqrt(S) * P * Sqrt(S)
    Pop = sum([M[i, i] for i in Fragment.FragmentBasisIterator()])

    assert abs(Pop.imag) < 1e-9
    return Fragment.Z - Pop.real


def LowdinFragmentDensity(P, S, Fragment):
    """Calculates Lowdin density fragment"""
    PI = Fragment.ProjectionMatrix()
    M = PI * Sqrt(S) * P * InvSqrt(S)
    return 0.5 * PurifyReal(M + M.H)



def MullikenFragmentDensity(P, Fragment):
    """Calculates Mulliken equivalent of the fragment density"""

    PI = Fragment.ProjectionMatrix()
    return 0.5 * (P * PI + PI * P)



def MullikenCharge(P, S, Fragment):
    """Calculates Mulliken charge for the fragment"""
    M = P * S
    Pop = sum([M[i, i] for i in Fragment.FragmentBasisIterator()])
    return Fragment.Z - Pop



