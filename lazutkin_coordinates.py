from mpmath import *
from domain import *

def C(domain):
    """ Returns the C coefficient for the Lazutkin paramterization.
    """

    integrand = lambda t: fmul(power(domain.radius(t), fdiv(-2, 3)), 
                               norm(domain.polar_gradient(t)))
    return power(quad(integrand, [0, 1]), -1)


def lazutkin_param_non_arc(domain, s):
    """ Returns the Lazutkin parametrization for domain evaluated at s.
        Generally s is assumed to be the parameterization by arc length.
    """

    integrand = lambda t: fmul(power(domain.radius(t), fdiv(-2, 3)), 
                               norm(domain.polar_gradient(t)))
    return fmul(C(domain), quad(integrand, [0, s]))

def lazutkin_weight(domain, x):
    """ Returns the Lazutkin weight of domain at the point x.
    """
    return power(fprod([2, C(domain), power(domain.radius(x), fdiv(1, 3))]), -1)