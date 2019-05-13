from mpmath import *
from domain import *

def lazutkin_param(domain, s):
    """ Returns the Lazutkin parametrization for domain evaluated at s.
        Generally s is assumed to be the parameterization by arc length.
    """

    integrand = lambda s_prime: power(domain.radius(normalize_coords(s_prime)), 
        fdiv(-2, 3))
    C = power(quad(integrand, [0, 1]), -1)
    return fmul(C, quad(integrand, [0, s]))

def lazutkin_weight(domain, x):
    """ Returns the Lazutkin weight of domain at the point x.
    """

    integrand = lambda s_prime: power(domain.radius(normalize_coords(s_prime)), 
        fdiv(-2, 3))
    C = power(quad(integrand, [0, 1]), -1)
    return power(fprod([2, C, power(domain.radius(x), fdiv(1, 2))]) , -1)