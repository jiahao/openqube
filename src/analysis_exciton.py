#!/usr/bin/env python

"""
Parses Q-Chem output files with FED calculations

For FED calculations, we grep TDDFT calculations
for the energies and transition dipoles of bright
excited states
"""

from numpy import dot
from numpy.linalg import norm

SpinMultiplicity = {
    "Singlet":0,
    "Triplet":2,
}

def nearest(x, things):
    "Returns the element in an iterable things closest to the querying x"
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

