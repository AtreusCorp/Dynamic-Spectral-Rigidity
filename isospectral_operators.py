from mpmath import *
from domain import *
from billiards import *
from lazutkin_coordinates import *

def fourier_basis(j):
    """ Returns the function e_j(x) = cos (2 pi j x) corresponding to the 
        fourier basis.
    """

    return lambda x: cos(fprod([2, pi, j, x]))

def l_q(domain, q, function):
    """ Returns the result of evaluating the functional l_q. (TODO: SAY THIS BETTER).
    """
    
    if (q == 0):
        integrand = lambda x: fprod([fdiv(function(x), domain.radius(x)), 
                                     norm(domain.polar_gradient(x))])

        # Here we make the assumption that the length of the boundary of
        # domain is 1
        return quad(integrand, [0, 1])

    elif (q == 1):
        return function(0)

    orbit = compute_q_bounce_path(domain, q)
    summands = [fmul(function(point[0]), sin(point[1])) for point in orbit]
    return fsum(summands)

def l_q_lazutkin(domain, q, function):
    """ Returns the output of the linearized map modified by the Lazutkin 
        weight.
    """

    arg_pre_scaled = lambda x: fdiv(function(x), lazutkin_weight(domain, x))
    scaled_function = lambda x: arg_pre_scaled(lazutkin_param_non_arc(domain, x))

    return l_q(domain, q, scaled_function)

def T(domain, function, precision):
    """ Returns a list of length precision given by the linearized isospectral 
        operator.
    """

    output = []
    q = 0

    while (q < precision):
        output.append(l_q(domain, q, function))
        q += 1
    return output

def T_lazutkin(domain, function, precision):
    """ Returns a list of length precision given by the linearized isospectral 
        operator modified by the Lazutkin weight.
    """

    output = []
    q = 0

    while (q < precision):
        output.append(l_q_lazutkin(domain, q, function))
        q += 1
    return output

def T_lazutkin_matrix(domain, i, j):
    """ Returns the ith element of T_lazutkin(e_j).
    """

    return l_q_lazutkin(domain, i, fourier_basis(j))

def P_star(function):
    """ Returns the function obtained by subtracting the integral from 0 to 1
        of function from function.
    """

    integral = quad(function, [0, 1])
    return lambda x: function(x) - integral

def b_bullet(precision):
    """ Compute the tuple (0, 0, 1/4, ... , 1 / q^2, ...) up to at least 
        precision dimensions.
    """

    b = [0, 0]
    q = 2

    while (q < precision):
        b.append(fdiv(1, power(q, 2)))
        q += 1
    return b

def operator_norm(matrix, gamma, max_j, max_q):
    """ Returns the gamma operator norm of matrix, summing up to max_j and
        considering the sup up to max_q. Assumed that matrix is a function 
        accepting two arguments i,j and not an array () for efficiency.
    """

    max_j_sum = -1
    q = 0

    while(q < max_q):
        temp_j_sum = nsum(lambda j: fprod([power(q, gamma), power(j, -gamma),
                                           matrix(q, j)]), [1, max_j])
        max_j_sum = temp_j_sum if temp_j_sum > max_j_sum else max_j_sum
        q += 1
    return max_j_sum