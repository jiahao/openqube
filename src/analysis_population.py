#!/usr/bin/env python
#######################################
# Functions for population analysis

from scipy import asmatrix, matrix, diag, real, sqrt, finfo
from numpy.linalg import eig, inv, norm

_epsilon = sqrt(finfo(float).eps)

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
    return 0.5 * real(M + M.H)

def MullikenFragmentDensity(P, Fragment):
    """Calculates Mulliken equivalent of the fragment density"""
    PI = Fragment.ProjectionMatrix()
    return 0.5 * (P * PI + PI * P)

def MullikenCharge(P, S, Fragment):
    """Calculates Mulliken charge for the fragment"""
    M = P * S
    Pop = sum([M[i, i] for i in Fragment.FragmentBasisIterator()])
    return Fragment.Z - Pop

def isHermitian(M, tol):
    MM = asmatrix(M)
    return norm(MM - MM.H) < tol

def Pow(M, power = 1):
    """Calculates the principal power of a square Hermitian matrix"""
    assert isHermitian(M, tol = 1e-14)
    D, U = eig(M)
    U = matrix(U)
    E = [x ** power for x in D]
    E = matrix(diag(E))
    Mpow = U * E * inv(U)
    return Mpow

def Sqrt(M):
    """Calculates the principal square root of a square Hermitian matrix"""
    return Pow(M, 0.5)

def InvSqrt(M):
    """Calculates the principal inverse square root of a square Hermitian matrix"""
    return Pow(M, -0.5)

def Inv(M):
    """Calculates the principal inverse of a square matrix"""
    return Pow(M, -1)

