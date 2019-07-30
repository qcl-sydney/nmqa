'''
Implements Padua interpolant from:

[1] Caliari M, De Marchi S, Sommariva A, Vianello M. Padua2DM: fast interpolation and
cubature at the Padua points in Matlab/Octave. Numerical Algorithms. 2011 Jan 1;56(1):45-60.

'''
import numpy as np

from pdpoints import dims_padua_set, T_n, padua_cubature_weights

def Tn_(n, x_i):
    '''Return the  scaled Chebyshev polynomial of degree n'''

    return T_n(n, x_i) * np.sqrt(2)


def reproducing_kernel_v2(n, x, y):
    '''Return the reproducign kernel K_n(x,y) using Eq. (3) in [1]. '''

    ans = 0.

    for idx_k in range(n + 1):

        for idx_j in range(idx_k + 1):

            ans += Tn_(idx_j, x[0]) * Tn_(idx_k - idx_j, x[1]) * Tn_(idx_j, y[0]) * Tn_(idx_k - idx_j, y[1])

    return ans


def fundamental_L_B_v2(n, X, B):
    ''' Return L_B(X), the coefficient for a fundamental Lagrange polynomial
    interpolant evalulated at padua point B for arbitrary point X via Eq.2 in [1]
     '''

    if n == 0:
        weight_b = 2.
    if n > 0:
        weight_b = padua_cubature_weights(n, B)

    kernel = reproducing_kernel_v2(n, X, B) - (0.5 * Tn_(n, B[0]) * Tn_(n, X[0]))
    return weight_b * kernel


def f_interpolant_v2(n, x1, x2, padua_points, data_points):
    ''' Interpolant using fundamental lagrange polynomials at the Padua points
    Reference:
        BIVARIATE LAGRANGE INTERPOLATION AT THE PADUA POINTS: THE IDEAL THEORY APPROACH (2007)
    '''
    L_B_vector = np.asarray([ fundamental_L_B_v2(n, [x1, x2], B) for B in padua_points])
    f_interpolant = np.sum(L_B_vector * data_points)    


    return f_interpolant
