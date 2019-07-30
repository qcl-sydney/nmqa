'''
Implements Padua interpolant from:

[1] Caliari M, De Marchi S, Sommariva A, Vianello M. Padua2DM: fast interpolation and
cubature at the Padua points in Matlab/Octave. Numerical Algorithms. 2011 Jan 1;56(1):45-60.

'''
import numpy as np

from pdpoints import padua_cubature_weights, T_n
from pdweights import Tn_

def chebyshevgausslobatto(n):
    '''
    Return  Chebyshev Gauss Lobatto points, dims: n + 1
    in p. 2 of [1].
    '''
    z_n_j = np.cos(np.arange(0., n + 1.)*np.pi / float(n))
    return z_n_j

def calc_t_matrix(n, c_vector):
    '''
    Return the rectangular Chebyshev matrix where
        n is the order of Padua points
        c_vector is a vector of chebyshevgausslobatto points
    in Eq. 6 of [1]

    '''
    # theory paper [1]
    # t_matrix = np.row_stack((Tn_(idx_n, c_vector) for idx_n in range(n+1)))

    # Cheb2fun
    t_matrix = np.row_stack((T_n(idx_n, c_vector) for idx_n in range(n+1)))

    return t_matrix

def calc_padua_cgl(order):

    ''' Return Padua points, Padua weights ad the  Chebyshev Gauss Lobatto points.
    # Points adhere to a consistent order of how the two sets for n+1, n+2 are
    # merged to from Chebyshev-like grid.
    '''
    # # Make a mask to pick out Padua points
    # g_idx = np.zeros((order+1)*(order+2), dtype='int')
    # g_idx[::2] = 1

    # # Mask differs by odd and even n
    # if np.mod(order, 2) == 1: # odd
    #     g_map = g_idx.reshape(order + 1, order + 2)
    # if np.mod(order, 2) == 0: # even
    #     g_map = g_idx.reshape(order + 2, order + 1).T

    # # Generate chebyshev gauss lobatto point vectors
    # z_s_vector = chebyshevgausslobatto(order + 1)
    # z_r_vector = chebyshevgausslobatto(order)
    # zsgrid, zrgrid = np.meshgrid(z_s_vector, z_r_vector)

    # # Pick out Padua points
    # pd_x1 = -1.0 * zrgrid[g_map == 1]
    # pd_x2 = -1.0 * zsgrid[g_map == 1]

    # CHEB2FUN
    g_idx = np.ones((order+1)*(order+2), dtype='int')
    g_idx[::2] = 0

    # Mask differs by odd and even n
    if np.mod(order, 2) == 1: # odd
        g_map = g_idx.reshape(order + 1, order + 2)
    if np.mod(order, 2) == 0: # even
        g_map = g_idx.reshape(order + 2, order + 1).T

    # Generate chebyshev gauss lobatto point vectors
    z_s_vector = chebyshevgausslobatto(order + 1)[::-1]
    z_r_vector = chebyshevgausslobatto(order)[::-1]
    zsgrid, zrgrid = np.meshgrid(z_s_vector, z_r_vector)

    # Pick out Padua points
    pd_x1 = -1.0 * zrgrid[g_map == 1]
    pd_x2 = -1.0 * zsgrid[g_map == 1]

    padua_points = zip(pd_x1, pd_x2)

    weights_cgl = np.asarray([padua_cubature_weights(order, item) for item in padua_points])

    return padua_points, weights_cgl, g_map, z_r_vector, z_s_vector

def calc_g_matrix(n, weights_cgl, g_map, f_data):
    '''
    Return matrix computed corresponding to the Chebyshev-like grid, p. 4 of [1].
    '''
    g_matrix = np.zeros(((n + 1), (n + 2)))

    # theory paper [1]
    # g_matrix[g_map == 1] = weights_cgl * f_data

    # Cheb2fun
    g_matrix[g_map == 1] = weights_cgl * f_data * 4.0

    return g_matrix

def pd_coefficients_matrix(n, f_data):
    '''
    Return $C_0$ coefficient matrix from Eq. 10  of [1] using MM approach.

    '''

    points, weights, g_map, c_vector_r, c_vector_s = calc_padua_cgl(n)

    g_matrix = calc_g_matrix(n, weights, g_map, f_data)

    # theory paper [1]
    # t_matrix_r = calc_t_matrix(n, c_vector_r)
    # t_matrix_s = calc_t_matrix(n, c_vector_s)
    # c_matrix = np.linalg.multi_dot([t_matrix_r, g_matrix, t_matrix_s.T])

    #cheb2fun
    t_matrix_r = np.cos(np.outer(range(n+1), range(n+1))*np.pi/float(n))
    t_matrix_s = np.cos(np.outer(range(n+2), range(n+2))*np.pi/float(n+1))
    c_matrix = np.linalg.multi_dot([t_matrix_s, g_matrix.T, t_matrix_r])
    c_matrix[0, :] = 0.5*c_matrix[0, :]
    c_matrix[:, 0] = 0.5*c_matrix[:, 0]
    c_matrix[0, -1] = 0.5*c_matrix[0, -1]

    coefficient_matrix = np.fliplr(np.triu(np.fliplr(c_matrix))) #  Eq. 10 # extra flipr? No.
    # theory paper [1]
    # coefficient_matrix[-1, 0] = coefficient_matrix[-1, 0] * 0.5  #  Eq. 10
    # # applies here. Agrees with author MM code

    return coefficient_matrix[0:n+1, :]


def pd_interpolant(n, X_datapts, X_testpts):
    '''
    Return Lnf(X); according to Eq. 11 and 12 of [1].
    '''

    # n +1 by n + 1 matrix of coefficients incorporating data
    c_0_matrix = pd_coefficients_matrix(n, X_datapts)

    t_x0_matrix = calc_t_matrix(n, X_testpts[0].flatten())  # vector of x coordinates of test points
    t_x1_matrix = calc_t_matrix(n, X_testpts[1].flatten()) # vecotr of y coordinates of test points

    # theory paper [1]
    # interpolant = np.linalg.multi_dot([t_x0_matrix.T, c_0_matrix, t_x1_matrix])

    # cheb2fun
    interpolant = np.linalg.multi_dot([t_x1_matrix.T, c_0_matrix, t_x0_matrix])

    return interpolant


