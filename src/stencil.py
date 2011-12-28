#!/usr/bin/env python
"""@package Stencil

Stencils for finite difference
Created on Jan 20, 2011

@author: Jiahao Chen
"""

import scipy as sc
import scipy.linalg as la
import logging as log

class Stencil:
    """
    Data structure containing information on finite derivative coefficients

    The finite difference approximation for a derivative can be written in the
    general stencil form

    @f[
    f^{(n)}(x) \approx \sum_{k \in S} c_k f(x + k \delta x)
    @f]

    This class stores information about the stencil set @f$ S @f$ and the
    coefficient set @f$ {c_k} @f$, and the order of differentiation @f$ n @f$.
    """
    def __init__(self, Coefficients = ((0, 0.0),), order = 1, Base = 1):
        """
        @param Coefficients An iterable of (stencil_point, unnormalized
        weight) tuples. Stencil_point is an integer and the weight
        is a float. Default: ((0, 0.0),)
        @param order Order of derivative (first = 1, second = 2, ...)
        @param Base Normalization factor for the weights. If given, will divide
               all the weights. Default: 1.0
        """

        if Base == 1:
            self.coeffs = list(Coefficients)
        elif Base != 0:
            ## An iterable of (stencil_point, weight) tuples.
            ## stencil_point is an integer and the weight is a float.
            ## is a float.
            ## The collection of stencil_points defines the stencil set
            ## @f$ S @f$.
            self.coeffs = [(n, c / float(Base)) for (n, c) in Coefficients]
        else: #Base == 0
            raise ValueError, 'Divide by zero in initialization of stencil.'

        ## Order of derivative (first = 1, second = 2, ...)
        self.order = order



    def AddZero(self):
        """
        Adds the origin point with zero weight if absent.
        """
        if 0 not in [n for (n, _) in self.coeffs]:
            self.coeffs.append((0, 0))


    def GetIndex(self, nn):
        """
        Find the internal index in which nn is located.

        @param nn Element of stencil set to locate

        @returns Internal index. If not found, returns @c None.

        @todo This is inelegant
        """
        idx = [i for i, (n, _) in enumerate(self.coeffs) if n == nn]
        if len(idx) == 0:
            return None
        elif len(idx) > 1:
            assert False, 'Error, more than one found'
            return idx
        else:
            return idx[0]


    def GetNumPoints(self):
        """
        Computes number of points in Stencil.
        """
        return len(self.coeffs)


    def ApplyToFunction(self, f, x, h, *args):
        """
        Applies Stencil to the scalar function f at point x with grid spacing h
        
        @param f Function to be called, e.g. with f(x, *args)
        @param x Origin point to evaluate f at
        @param h Size of stencil grid
        @param args Optional additional arguments to f

        @returns Numerical derivative
        
        @todo Currently only does this in one dimension
        """

        x = sc.asarray(x)
        h = sc.asarray(h)
        f_x = f(x, *args)
        Df = sc.zeros(f_x.shape + x.shape)
        log.info('Generating numerical derivative of shape ' +str(Df.shape))
        #import itertools
        if h.ndim == 0:
            # itertools to make multiindices out of cartesian products of
            # 1-d ranges
            for n, c in self.coeffs:
                for x_multiindex in range(x.shape[0]): 
                    tmp = x[x_multiindex]
                    x[x_multiindex] += n * h
                    Df[:, x_multiindex] += c * f(x, *args)
                    x[x_multiindex] = tmp
                #for Df_multiindex in itertools.product(\
                #                 *itertools.imap(range, Df.shape)):
                #    log.error("Stencil point %d: Weight: %f " % (n, c) + \
                #              "Function value: " + str(f_x[f_multiindex]))
                #    log.error("Current multiindex: " + str(Df_multiindex))
                #    f_multiindex = Df_multiindex[:f_x.ndim]
                #    x_multiindex = Df_multiindex[f_x.ndim:]
                #    Df[Df_multiindex] += c * f_x[f_multiindex]
        else:
            #@todo extend this to other possibilities
            assert False, 'Not implemented, dim(h) = %d' % h.ndim

        return Df * la.norm(h) ** -self.order



    def CheckDeriv(self, f, df, x, h, normthresh = None, *args):
        """
        Checks the derivative df using this Stencil.

        @param f Function to check. Takes as arguments (x, *args)
        @param df The function to be tested as a derivative of f,
               which takes as arguments (x, *args),
               or the value of the derivative function evaluated at this point

        @param x Origin point to evaluate f and df at
        @param h Size of stencil grid
              
        @param normthresh Threshold of a norm to accept. Default: None
        
        @param args Any additional arguments required by f and df.

        @returns Norm of the difference if normthresh is @c None, or
                 the entire difference vector if normthresh is 0, or
                 @c True if the norm of the difference does not exceed
                 normthresh, otherwise @c False
        """

        try:
            dDf = df - self.ApplyToFunction(f, x, h, *args)
        except TypeError:
            dDf = df(x, *args) - self.ApplyToFunction(f, x, h, *args)

        if normthresh == 0:
            return dDf

        normdDf = la.norm(dDf)

        if normthresh == None:
            return normdDf
        elif normthresh > 0:
            return normdDf < normthresh



    def ApplyToFunctionOnGrid(self, f_x, h):
        """
        Applies Stencil to a function already evaluated on the necessary
        grid points in f_x with grid spacing h

        @param f_x Grid of evaluated function values

        @param h Size of stencil grid

        @returns Numerical derivative
        
        @todo Currently only does this in one dimension
        """

        Df = sum([c * f_x[i] for i, (_, c) in enumerate(self.coeffs)])

        return Df * h ** -self.order



####################
# Two-point stencils
class FirstOrderForwardDifferenceStencil(Stencil):
    """First-order classical forward difference Stencil"""
    def __init__(self):
        """
        Constructor.
        """
        Stencil.__init__(self,
                         Coefficients = ((0, -1), (1, 1)),
                         order = 1)



######################
# Three-point stencils

class FirstOrderCentralDifferenceStencil(Stencil):
    """First-order classical central difference Stencil"""
    def __init__(self):
        """
        Constructor.
        """
        Stencil.__init__(self,
                         Coefficients = ((-1, -1), (1, 1)),
                         Base = 2,
                         order = 1)



#####################
# Five-point stencils

class FirstOrderFivePointCentralDifferenceStencil(Stencil):
    """
    First-order five-point classical central difference Stencil
    """
    def __init__(self):
        """
        Constructor.
        """
        Stencil.__init__(self,
                         Coefficients = ((-2, 1), (-1, -8), (1, 8), (2, -1)),
                         Base = 12,
                         order = 1)



class FirstOrderFivePointLanczosDifferenceStencil(Stencil):
    """
    First-order five-point low-noise Lanczos difference Stencil
    """
    def __init__(self):
        """
        Constructor.
        """
        Stencil.__init__(self,
                         Coefficients = ((-2, -2), (-1, -1), (1, 1), (2, 2)),
                         Base = 10,
                         order = 1)



class FirstOrderFivePointHoloborodkoDifferenceStencil(Stencil):
    """
    First-order five-point smooth Holoborodko difference
    @see Pavel H., http://www.holoborodko.com/pavel/numerical-methods/numerical\
-derivative/smooth-low-noise-differentiators

    """

    def __init__(self):
        """
        Constructor.
        """
        Stencil.__init__(self,
                         Coefficients = ((-2, -1), (-1, -2), (1, 2), (2, 1)),
                         Base = 8,
                         order = 1,
                         )


class FirstOrderSevenPointDifferenceStencil(Stencil):
    """Chiao's magic first-order seven-point difference Stencil
    ripped from foptutils9.f1d7p()
    """
    def __init__(self):
        """
        Constructor.
        """
        Stencil.__init__(self,
                         Coefficients = ((3, 1), (2, -9), (1, 45), (-1, -45),
                                         (-2, 9), (-3, -1)),
                         Base = 60,
                         order = 1)

