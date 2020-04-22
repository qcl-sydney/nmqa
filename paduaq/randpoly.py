'''
Implements a scaled and shifted random bivariate polynomial of degree $n$ such that values
over relevant points sets are [0, \pi]

References:

[1] Alves Diniz, Marcio and Salasar, Luis Ernesto and Bassi Stern, Rafael. `Positive Polynomials on closed boxes`, arXiv preprint arXiv:1610.01437 (2016). 

[2] Gasca, Mariano, and Thomas Sauer. `Polynomial interpolation in several variables.` Advances in Computational Mathematics 12.4 (2000): 377.

'''

import numpy as np, sys
from scipy.special import binom


NEARESTPADUA={}
NEARESTPADUA["regfine"] = 9 #10 # 9x9 grid but 54 qubits
NEARESTPADUA["regcoarse"] = 4 # 4x4 grid with 16 qubits
NEARESTPADUA["reg4"] = 1 # 2x2 grid with 4 qubits
NEARESTPADUA["reg9"] = 3 # 3x3 grid with 9 qubits
NEARESTPADUA["reg36"] = 7 #5 # 6x6 grid with 36 qubits
NEARESTPADUA["reg25"] = 6 #5 # 5x5 grid with 25 qubits (shrunk)

def dims_coeffcients(n, m):
    '''Return number of coefficients for polynomial of degree n and
    inderminate variables m.'''
    
    return int(binom(n + m, n ))


def random_coeffients(n, m, mval=0.5):
    '''Return a set of random variates uniformly distrbuted on 
       the real half open interval [-mval, mval).
       Number of random variates is determined by n, m.
    '''
    dims =  dims_coeffcients(n, m)
    coeff = np.random.uniform(low=-mval, high=mval, size=dims)
    return coeff
    
    
def bivariate_lattice_nodes(n):
    ''' Return lattice node point set for m=2  bivariate polynomials 
    of degree n'''
    
    X, Y = np.meshgrid(np.arange(n+1),np.arange(n+1))
    alpha_1 = np.fliplr(np.triu(np.fliplr(X))).flatten()
    alpha_2 = np.fliplr(np.triu(np.fliplr(Y))).flatten()
    alpha_set = set(zip(alpha_1,alpha_2))

    for point in alpha_set:
        if point[0] + point[1] > n:
            print "Error", point
        if point[0] + point[1] <= n:
            pass
    
    if dims_coeffcients(n, 2) != len(alpha_set): 
        print "Error: dims of alpha point set does not match expected value "
        raise RuntimeError

    return alpha_set

def bivariate_polynomial(x1, x2, n, a_set):
    '''Return the value of a bivariate polynomial on 
    the point x=(x1, x2).
    Bivarate polynomial is of degree n and coefficients a_set '''
    
    alphaset = list(bivariate_lattice_nodes(n))
    dims = dims_coeffcients(n, 2)

    output=0
    for idx_pt in range(dims):
        alpha_1, alpha_2 = alphaset[idx_pt]
        a_coeff = a_set[idx_pt]
        output += a_coeff * (x1**alpha_1) * (x2**alpha_2)
        
    return output


def evaluate_bivariate_poly(points, n, a_set):
    '''Return values of a bivariate polynomial of degree n 
        and coefficients a_set on some list of points.'''
    
    f_set=[]
    p_set=[]
    for point in points:
        x1, x2 = point
        fval = bivariate_polynomial(x1, x2, n, a_set)
        f_set.append(fval)
        p_set.append(point)
    
    return p_set, np.asarray(f_set)
    
    
def scale_shift_polynomial(points, n, a_set):
    '''Return scale and shift parameters so that output of a
    bivariate polynomial is [0, \pi ] on some list of points. '''
    
    fvals = evaluate_bivariate_poly(points, n, a_set)[1]
    delta = max(fvals) - min(fvals)
    scale = np.pi / delta
    shift = - min(fvals)
    
    return scale, shift

def get_scaled_random_bivariatepoly(points, n, a_set):
    '''Return values of a scaled and shifted polynomial of degree $n$
    and real coefficients a_set such that all values are
    within [0, \pi] over the list of points'''
    
    scale, shift = scale_shift_polynomial(points, n, a_set)
    sfvals = scale*(evaluate_bivariate_poly(points, n, a_set)[1] + shift)
    
    return sfvals
