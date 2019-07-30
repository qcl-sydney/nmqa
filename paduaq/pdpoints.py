import numpy as np

# ===============================================================================================
# Global Padua Properties
# ===============================================================================================
def generating_curve(n, t):
    '''Return the x1 and x2 coordinates of the point x on a Padua generating curve
        of order n and angle t
    Input:
        n - max degree of polynomial in space of polynomials of degree at most n in two variable
        t - angle; or spacing length along the generating curve.
            Calculated by calc_padua_points() as equi-spaced intervals to make Padua points

    Reference:
        Bos et. al. Bivariate Lagrange interpolation at the Padua points:
            the generating curve approach. (2006) - incorrect
        Padua2DM: fast interpolation and cubature at the Padua
            points in Matlab/Octave - CORRECT
        '''

    x2 = -1.0 * np.cos(n*t)
    x1 = -1.0 * np.cos((n + 1.)*t)
    return np.asarray([x1, x2])

def dims_padua_set(n):
    '''Return the dimensions of the Padua set of order n.
    Input:
        n - max degree of polynomial in space of polynomials of degree at most n in two variable.
    Reference:
        Bos et. al. Bivariate Lagrange interpolation at the Padua points:
        the generating curve approach. (2006)
    '''

    return (n + 2.)*(n + 1.) / 2.

# ===============================================================================================
# Padua Point Set Generation  - Ideal Theory Approach (modified)
# Reference: BIVARIATE LAGRANGE INTERPOLATION AT THE PADUA POINTS: THE IDEAL THEORY APPROACH (2007)
# ===============================================================================================

def calc_padua_points_v2(n, sigfig=6):

    '''
    Return Pad_n, the set of Padua points on the square [-1,1]^2, and their weights
    Uses a modified Ideal Theory Approach.

    Input:
        n - max degree of polynomial in space of polynomials of degree at most n in two variable

    Reference:
        BIVARIATE LAGRANGE INTERPOLATION AT THE PADUA POINTS: THE IDEAL THEORY APPROACH (2007)
        De Marchi, Stefan. Padua points:  genesis, theory, computation and applications (2014)
        : http://www.math.iit.edu/~Meshfree-methods-seminar/presentations/talk_20140115_stefano.pdf
    '''
    dims = int(dims_padua_set(n))

    pointset = generate_padua_points(n, sigfig=sigfig)
    weights = [padua_cubature_weights(n, pnt) for pnt in pointset]

    return pointset, weights

def padua_cubature_weights(n, x):
    ''' Wrapper function for get_weight() in order to troubleshoot theory typos in Padua literature.
    '''

    global_weight = float(n * (n + 1.0))

    # return 1.0 / (get_weight(n, x) * global_weight) # IDEAL THEORY approach w = 1/K*(x, x) in Prop 3.3; K*(x, x) = get_weight(n, x) * global_weight
    return get_weight(n, x) / global_weight #  w_A from Generating Curve Approach In Th.1 ; equivalently w_epsilon in De Marchi talk (2014) p33


def get_weight(n, x, vertex_weight=0.5,
                  edge_weight=1.,
                  interior_weight=2.0):

    ''' Return the weight of Padua point x, by classifying it as a boundary,
    vertex or interior point.

    Notes:
        Vertex idenitifcation from Ideal Theory approach (above Prop 3.3)
            doesnt work, so we try an alternative technique.
        Ref: BIVARIATE LAGRANGE INTERPOLATION AT THE PADUA POINTS: THE IDEAL THEORY APPROACH (2007)
    '''

    # Vertex
    if abs(x[0]) == abs(x[1]): # modified approach to Ideal Theory approach
        return vertex_weight

    # Boundary
    for x_coord in [x[0], x[1]]:
        if abs(x_coord) == 1.0:
            return edge_weight

    # Interior
    return interior_weight


def generate_padua_points(n, sigfig=6):

    '''
    Return Pad_n, the set of Padua points on the unit square [-1,1]^2.

    Input:
        n - max degree of polynomial in space of polynomials of degree at most n in two variable

    Notes:
        Cross-checked with Transformed Generating Gurve approach (my code) and padua.py [authors]
        Eqns (1.1) and (1.2) in ideal theory approach - modified floor function to cieling
            function to compute j-values in (1.1).
        Ref: BIVARIATE LAGRANGE INTERPOLATION AT THE PADUA POINTS: THE IDEAL THEORY APPROACH (2007)
    '''

    # jvalues  = np.arange(1, np.floor(n * 0.5) + 1. + 1.)
    # ## Eqns (1.1, 1.2) in Ideal theory approach paper
    # ## Generates correct point set for even n. But odd n, missing points for some k values, with x2 = -1.
    # ## Number of missing values = np.ceil(n/2), n odd
    # ## Corresponds to left edge boundary points (x2=-1) not being computed for odd n

    jvalues = np.arange(1, np.ceil(n * 0.5) + 1. + 1.)  # modified Ideal Theory approach
    # Generates correct point set for odd n, and even n, but with duplicate points
    # ## Duplicates are removed by rounding to sigfig decimal points and using set()

    dims_j = len(jvalues)
    dims_k = n + 1
    dims_total = int(dims_padua_set(n))

    point_set = np.zeros((dims_total, 2))

    list_of_tuples = []
    for idx_k in np.arange(dims_k):  # 0 <= k <= n

        x1_values = np.cos(idx_k * np.pi / n) * np.ones(dims_j)

        if idx_k % 2 != 0: # ( k is odd)
            x2_values = np.cos((2 * jvalues - 2) * np.pi / (n + 1.)) # broadcasting

        if idx_k % 2 == 0: # ( k is even)
            x2_values = np.cos((2 * jvalues - 1) * np.pi / (n + 1.)) # broadcasting

        # modified Ideal Theory approach
        list_of_tuples += zip(np.round(x1_values, sigfig), np.round(x2_values, sigfig))

    # modified Ideal Theory approach
    remove_duplicates = list(set(list_of_tuples))
    point_set[:, 0], point_set[:, 1] = np.asarray(zip(*remove_duplicates))

    return point_set

# ===============================================================================================
# Padua Point Set Generation  - Generating Curve Approach (modified)
# Reference: Bos et. al. Bivariate Lagrange interpolation at the Padua points: the generating curve approach. (2006)
# ===============================================================================================

def padua_index_weights(n, 
                        edge_weight=1., 
                        vertex_weight=0.5,
                        interior_weight=2.):
    
    '''Return the pairs of indices (j,m) and their weight for each point in Padua set of order n
    
    Input:
        n - max degree of polynomial in space of polynomials of degree at most n in two variable
    Reference:
        Bos et. al. Bivariate Lagrange interpolation at the Padua points: the generating curve approach. (2006)
    
    '''
    
    dims = int(dims_padua_set(n))
    global_weight = (n * (n + 1.))
    
    # VERTEX
    vertex_0 = np.zeros((1, 3))
    vertex_0[:, 2] = vertex_weight * global_weight
    
    # VERTICAL EDGES
    vert_edge_m = None
    interior_indices = None
    if n > 1:
        vert_edge_m = np.zeros((n-1, 3))
        vert_edge_m[:, 1] = np.arange(1, n)
        vert_edge_m[:, 2] = edge_weight * global_weight
               
        interior=[]
        for i in range(n):
            j = i + 1.0
            m_max = n - j 

            if m_max > 0 :
                for item in range(1, int(m_max + 1)):
                    interior.append([j, float(item)])
        
        interior_indices = np.zeros((len(interior), 3))
        interior_indices[:,0:2] = np.asarray(interior)
        interior_indices[:, 2] = interior_weight * global_weight
    
    # VERTEX
    vertex_n = None
    
    # HORIZONTAL EDGES
    hort_edge_j= None
    if n > 0:
        vertex_n = np.zeros((1, 3))
        vertex_n[0,1] = n
        vertex_n[0,2] = vertex_weight * global_weight
        
        #print ("interior_weight: ", vertex_weight * global_weight)

        hort_edge_j = np.zeros((n, 3))
        hort_edge_j[:, 0] = np.arange(1, n + 1)
        hort_edge_j[:, 2] = edge_weight * global_weight
        
        #print ("interior_weight: ", edge_weight * global_weight)
        
        
    # STACK VERTICES, VERTICAL EDGES, HORIZONTAL EDGES, INTERIOR POINTS
    index_set = np.vstack([vertex_0])
    for item in [vertex_n, vert_edge_m, hort_edge_j, interior_indices]:
        if item is not None:
            index_set = np.vstack([index_set, item])

    return index_set

def calc_padua_points(n):

    '''
    Return Pad_n, the set of Padua points on the square [-1,1]^2, and their weights
    Uses a Generating Curve Approach.

    Input:
    n - max degree of polynomial in space of polynomials of degree at most n in two variable
    '''
    index_weights = padua_index_weights(n)
    pts = int(dims_padua_set(n))
    padua_points = np.zeros((pts, 2))

    if n > 0:

        for idx_pt in range(pts):

            j, m = index_weights[idx_pt, 0:2]
            arg = ((j * n) + m * (n + 1.)) * np.pi / (n * (n + 1.))
            padua_points[idx_pt, :] = generating_curve(n, arg)

    return padua_points, index_weights[:, 2]

def transform_points(padua_points):
    '''
    Modifications to Generating Curve Approach
    Empirical plots show that Padua points generated by Bos et. al. in 2006
    correspond to a global -1.0 (rotation) and swapped x1, x2 coordinates (reflection).

    '''

    transformed_points = np.zeros_like(padua_points)
    transformed_points[:, 0] = -1.0 * padua_points[:, 1]
    transformed_points[:, 1] = -1.0 * padua_points[:, 0]

    return transformed_points



# ===============================================================================================
# Fundamental Lagrange Polynomial Calculations on Padua Point Set (VIA IDEAL THEORY APPROACH)
#
# Interpolation functions below do not use Padua weights (e.g. via padua_cubature_weights(n, x))
# Instead Lagrangian basis function are computed directly as: K*(x,y) / K(x,x).
# 
# References:
#   Bos et. al. Bivariate Lagrange interpolation at the Padua points: the generating curve approach. (2006)
#   BIVARIATE LAGRANGE INTERPOLATION AT THE PADUA POINTS: THE IDEAL THEORY APPROACH (2007)
# ===============================================================================================


def K_star(n, x, y):
    
    '''
    Reference:
    Proposition 3.1. in Bos et. al. BIVARIATE LAGRANGE INTERPOLATION AT THE PADUA POINTS: THE IDEAL THEORY APPROACH (2006)

    '''
    ans = reproducing_kernel(n, x, y) - T_n(n, x[0])*T_n(n, y[0])
    
    return ans

def reproducing_kernel(n, A, B):
    ''' Return the reproducing kernel for the inner product defined on the space of polynomials on a square.
    
    Reference:
        Bos et. al. Bivariate Lagrange interpolation at the Padua points: the generating curve approach. (2006)
        In particular, refer Lemma 2 for the form of the reproducing kernel for any two points A, B on the square.
    '''
    
    theta1 = np.arccos(A[0])
    theta2 = np.arccos(A[1])
    phi1 = np.arccos(B[0])
    phi2 = np.arccos(B[1])
    
    ans=0.
    ans += D_operator(n, theta1 + phi1, theta2 + phi2)
    ans += D_operator(n, theta1 + phi1, theta2 - phi2)
    ans += D_operator(n, theta1 - phi1, theta2 + phi2)
    ans += D_operator(n, theta1 - phi1, theta2 - phi2)
    
    return ans

def D_operator(n, alpha, beta):
    
    '''Helper function for reproducing kernel function.
    
    Reference: Lemma 2 in
        Bos et. al. Bivariate Lagrange interpolation at the Padua points: the generating curve approach. (2004)
    '''
    
    numer = np.cos(alpha * (n + 0.5)) * np.cos(alpha *  0.5)
    numer += -1.0 * np.cos(beta * (n + 0.5)) * np.cos(beta *  0.5)
    
    denom = np.cos(alpha) - np.cos(beta)
    
    ans  = 0.5 * numer / denom
    
    #print("D operator", ans)
    
    if np.isnan(ans):
        
        if numer == 0. and denom == 0. :
            ans = 1.
            
            #print("D operator reset to", ans)
            
        if denom == 0. and numer != 0:
            print ("Ah fuck")
    
    return ans


def fundamental_L_B(n, X, B):
    ''' Return L_B(X), the coefficient for a fundamental Lagrange polynomial
    interpolant evalulated at padua point B for arbitrary point X. 
    
    The polynomials L_B are indeed the fundamental Lagrange polynomials,i.e., they satisfy
    L_B(A) = \delta_{A,B};  for A, B in Pad_n; and \Delta = 1 if A==B else 0
    
    References: 
        Theorem 2 in Bos et. al. Bivariate Lagrange interpolation at the Padua points: the generating curve approach. (2006)
        Theorem 3.2. in Bos et. al. BIVARIATE LAGRANGE INTERPOLATION AT THE PADUA POINTS: THE IDEAL THEORY APPROACH (2007)
        Page 33. in De Marchi talk (2017)
    
    '''
       
    scalar_coeff_2 = K_star(n, X, B) / K_star(n, B, B) # IDEAL THEORY
    # scalar_coeff_2 = K_star(n, B, X) / K_star(n, B, B) # IDEAL THEORY change order
    # scalar_coeff_2 = K_star(n, X, B) * weight_B # should be equivalent according to IDEAL THEORY, GEN CURVE, and De Marchi Talk (2014)
    
    return scalar_coeff_2

def T_n(n, x_i):
    '''Chebyshev polynomial of the first kind of order n. Confirmed to be identifcal to Scipy.
    Supports broadcasting.
    
    Reference:
    Bos et. al. IVARIATE LAGRANGE INTERPOLATION AT THE PADUA POINTS: THE IDEAL THEORY APPROACH (2006)
    '''
    
    theta = np.arccos(x_i)
    ans = np.cos(n * theta)
    
    return ans

def f_interpolant(n, x1, x2, padua_points, data_points):
    ''' Interpolant using fundamental lagrange polynomials at the Padua points 
    
    Reference:
        BIVARIATE LAGRANGE INTERPOLATION AT THE PADUA POINTS: THE IDEAL THEORY APPROACH (2007)
    '''
    
    
    L_B_vector = np.asarray([ fundamental_L_B(n, [x1, x2], B) for B in padua_points])
    f_interpolant = np.sum(L_B_vector * data_points)    
    
    return f_interpolant
