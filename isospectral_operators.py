from mpmath import *
from domain import *
from billiards import *
from lazutkin_coordinates import *

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

    arg_pre_scaled = lambda x: fmul(lazutkin_weight(domain, x), function(x))
    scaled_function = lambda x: arg_pre_scaled(lazutkin_param(domain, 
        arc_length_coords(domain, x)))

    return l_q(domain, q, scaled_function)

def T(domain, function, precision):
    """ Returns a list of length precision given by the linearized isospectral 
        operator.
    """

    output = []
    q = 0

    while (q < precision):
        output.append(l_q(domain, q, function))
    return output

def T_lazutkin(domain, function, precision):
    """ Returns a list of length precision given by the linearized isospectral 
        operator modified by the Lazutkin weight.
    """

    output = []
    q = 0

    while (q < precision):
        output.append(l_q_lazutkin(domain, q, function))
    return output