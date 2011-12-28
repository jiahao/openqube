#!/usr/bin/env python
"""
Helper routines for dealing with matrices

Jiahao Chen, 2009-12
"""

import scipy
import scipy.linalg as linalg

def PrintMatrix(Mraw, title = '', color = False, diag = False):
    """Prints a matrix
    optional arguments:
    title prints an optional title
    color turns on a simple color coding scheme using ANSI escape sequences
    diag adds an asterisk to the end of a diagonal element to mark it
    """

    try:
        numrows, numcols = Mraw.shape
        M = Mraw
    except ValueError:
        M = scipy.asmatrix(Mraw.view())
        numrows = numcols = int((Mraw.shape[0]) ** 0.5)
        assert numrows * numcols == Mraw.shape[0], 'Wrong dimensions: (%d,) -X-> (%d, %d)' % (Mraw.shape[0], numrows, numcols)
        M.shape = (numrows, numcols)

    print
    print title
    def ansicolortable(x):
        if x < -1:
            escape = '91'
        elif x < -0.1:
            escape = '95'
        elif x < -0.01:
            escape = '35'
        elif x < 0.01:
            escape = '0'
        elif x < 0.1:
            escape = '34'
        elif x < 1:
            escape = '36'
        else:
            escape = '96'
        return '\033[' + escape + 'm'


    for i in range(numrows):
        for j in range(numcols):
            x = M[i, j]
            if color:
                print ansicolortable(x),

            if diag and i == j:
                print "%7.3f*" % x,
            else:
                print "%8.4f" % x,

            if color:
                print '\033[0m',
        print



def isHermitian(M, tol):
    MM = scipy.asmatrix(M)
    return linalg.norm(MM - MM.H) < tol


def Pow(M, power = 1):
    """Calculates the principal power of a square Hermitian matrix"""
    assert isHermitian(M, tol = 1e-14)
    D, U = linalg.eig(M)
    U = scipy.matrix(U)

    #Purify small imaginary roundoff errors
    #try:
    #    if max([x.imag for x in D]) < 1e-10:
    #        E = [x.real ** power for x in D]
    #    else:
    #        E = [x ** power for x in D]
    #except SyntaxError: #Real, not complex
    #    E = [x ** power for x in D]
    E = [x ** power for x in D]

    E = scipy.matrix(scipy.diag(E))
    Mpow = U * E * linalg.inv(U)

    return Mpow



def PurifyReal(M, tol = 1e-3):
    #try:
    #    for i in range(M.shape[0]):
    #        for j in range(M.shape[1]):
    #            assert M[i, j].imag < tol, 'PurifyReal: M[%d,%d].imag = %f > tol = %f' % (i, j, M[i, j].imag, tol)
    #            M[i, j] = M[i, j].real
    #except SyntaxError: #No complex type
    #    pass
    return scipy.real(M)


def Sqrt(M):
    """Calculates the principal square root of a square Hermitian matrix"""
    return Pow(M, 0.5)


def InvSqrt(M):
    """Calculates the principal inverse square root of a square Hermitian matrix"""
    return Pow(M, -0.5)



def Inv(M):
    """Calculates the principal inverse of a square matrix"""
    return Pow(M, -1)





def SymmetricPartAsVector(M):
    Na, Nb = M.shape
    assert Na == Nb, 'Non-square matrix'

    MM = 0.5 * (M + M.T)
    V = scipy.zeros(Na * (Nb + 1) / 2)
    mu = 0
    for i in range(Na):
        for j in range(i, Nb):
            V[mu] = MM[i, j]
            mu += 1
    return V


def MatrixFromSymmetricPartAsVector(V):
    N = int((0.25 + 2 * len(V)) ** 0.5 - 0.5)
    #print len(V), '-->', N, 'x', N
    M = scipy.zeros((N, N))
    mu = 0
    for i in range(N):
        for j in range(i, N):
            M[i, j] = V[mu]
            M[j, i] = V[mu]
            mu += 1
    return M



def Lyapunov(A, B, threshold = 1.0e-9, Verbose = False, CheckSane = True):
    """Solves the Lyapunov problem AX + X A.T = B"""
    n = A.shape[0]
    if n != A.shape[1] or n != B.shape[0] or n != B.shape[1]:
        print "Dimension error"
        raise ValueError

    D, U = linalg.eig(A)
    U = scipy.matrix(U)
    V = scipy.matrix(linalg.inv(U))

    if Verbose:
        print "In Lyapunov"
        print "-----------"
        print
        print "Eigendecomposition of A is:"
        PrintMatrix(U)
        print
        PrintMatrix(scipy.diag(D))
        print

    if CheckSane: #Runs sanity checks
        #Checks that the eigendecomposition is correct
        DD = scipy.matrix(scipy.diag(D))
        AA = U * DD * V
        pf = linalg.basic.norm(AA - A)

        if pf > threshold:
            print "Eigendecomposition wrong"
            print "Frobenius norm is", pf
            raise ValueError

    BB = V * B * V.T

    if Verbose:
        print "RHS transforms into"
        PrintMatrix(BB)
        print

    X = scipy.matrix(scipy.zeros((n, n)))

    nr = 0
    for i in range(n):
        for j in range(n):
            y = D[i] + D[j]
            if abs(y) > threshold:
                X[i, j] = BB[i, j] / y
            else:
                if Verbose:
                    print "X[%d, %d] was regularized - denominator =" % (i + 1, j + 1), y
                X[i, j] = 0
                nr += 1

    if Verbose:
        print "Transformed solution is"
        PrintMatrix(X)
        print

    if nr > 0:
        print "Warning: %d matrix elements fell below regularization threshold of %f" % (nr, threshold)

    X = U * X * U.T

    if CheckSane:
        Z = A * X + X * A.T - B
        zz = linalg.basic.norm(Z)
        if zz > threshold:
            print
            print "ERROR IN LYUPANOV SOLVER"
            print "------------------------"
            print
            print "Incorrect solution was computed"
            print
            #PrintMatrix(X)
            print
            print "Frobenius norm of error is ", zz
            print
            #raise ValueError

    return X



