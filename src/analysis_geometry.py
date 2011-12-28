#!/usr/bin/env python

"""
Parses Q-Chem output files with FED calculations

For FED calculations, we grep TDDFT calculations
for the energies and transition dipoles of bright
excited states
"""

from numpy import zeros

def MomentOfInertiaTensor(R, Weights = None, Center = True):
    "Moment of inertia tensor for a bunch of points, optionally weighted"
    I = zeros((3, 3))
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
                I[i, j] += weight * (p[i]-C[i]) * (p[j]-C[j])
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

