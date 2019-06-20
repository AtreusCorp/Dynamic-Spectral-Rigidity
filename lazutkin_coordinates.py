from mpmath import *
from domain import *

def C(domain):
    """ Returns the C coefficient for the Lazutkin paramterization.
    """

    integrand = lambda t: fmul(power(domain.radius_of_curv(t), fdiv(-2, 3)), 
                               norm(domain.polar_gradient(t)))
    return power(quad(integrand, [0, 1]), -1)


def lazutkin_param_non_arc(domain, s):
    """ Returns the Lazutkin parametrization for domain evaluated at s.
    """

    integrand = lambda t: fmul(power(domain.radius_of_curv(t), fdiv(-2, 3)), 
                               norm(domain.polar_gradient(t)))
    return fmul(C(domain), quad(integrand, [0, s]))

def lazutkin_param_non_arc_deriv(domain, s):
    """ Returns the Lazutkin parametrization for domain evaluated at s.
    """

    result = fmul(power(domain.radius_of_curv(t), fdiv(-2, 3)), 
                  norm(domain.polar_gradient(t)))
    return fmul(C(domain), result)

def lazutkin_param_non_arc_second_deriv(domain, s):
    """ Returns the Lazutkin parametrization for domain evaluated at s.
    """
    
    radius_of_curv = domain.radius_of_curv(s)
    radius_of_curv_deriv = domain.radius_of_curv_deriv(t)
    polar_gradient = domain.polar_gradient(t)
    polar_gradient_norm_deriv = domain.polar_gradient_norm_deriv(t)

    result = fadd(fprod([fdiv(-2, 3), power(radius_of_curv, fdiv(-5, 3)), 
                         radius_of_curv_deriv, norm(polar_gradient)]),
                  fmul(polar_gradient_norm_deriv, 
                       power(radius_of_curv, fdiv(-2, 3))))
    return fmul(C(domain), fnc(s))

def lazutkin_weight(domain, x):
    """ Returns the Lazutkin weight of domain at the point x.
    """
    
    return power(fprod([2, C(domain), 
                        power(domain.radius_of_curv(x), fdiv(1, 3))]), -1)