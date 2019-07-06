from mpmath import *
from domain import *

def C(domain):
    """ Returns the C coefficient for the Lazutkin paramterization.
    """

    integrand = lambda t: fmul(power(domain.radius_of_curv(t), fdiv(-2, 3)),
                               norm(domain.polar_gradient(t)))
    return power(quad(integrand, [0, 1]), -1)


def lazutkin_param_non_arc(domain, t):
    """ Returns the Lazutkin parameter for domain corresponding to the
        polar angle t.

    :param domain: the domain
    :param t: polar angle of the point under consideration

    """

    integrand = lambda s: fmul(power(domain.radius_of_curv(s), fdiv(-2, 3)),
                               norm(domain.polar_gradient(s)))
    return fmul(C(domain), quad(integrand, [0, t]))

def inv_lazutkin_param_non_arc(domain, x):
    """ Returns the inverse of the transformation taking x (in [0, 2 pi]) to the 
        corresponding point in the parameterization by arc length.
    """

    return findroot(lambda t: lazutkin_param_non_arc(domain, t) - x, 0.5)

def lazutkin_param_non_arc_deriv(domain, t):
    """Returns the derivative of the Lazutkin parameter for domain
    evaluated at t.

    FIXME: derivative with respect to what?

    :param domain: the domain
    :param t: polar angle of the point under consideration

    """

    result = fmul(power(domain.radius_of_curv(t), fdiv(-2, 3)),
                  norm(domain.polar_gradient(t)))
    return fmul(C(domain), result)

def lazutkin_param_non_arc_second_deriv(domain, t):
    """ Returns the Lazutkin parametrization for domain evaluated at t.
    """

    radius_of_curv = domain.radius_of_curv(t)
    radius_of_curv_deriv = domain.radius_of_curv_deriv(t)
    polar_gradient = domain.polar_gradient(t)
    polar_gradient_norm_deriv = domain.polar_gradient_norm_deriv(t)

    result = fadd(fprod([fdiv(-2, 3), power(radius_of_curv, fdiv(-5, 3)),
                         radius_of_curv_deriv, norm(polar_gradient)]),
                  fmul(polar_gradient_norm_deriv,
                       power(radius_of_curv, fdiv(-2, 3))))
    return fmul(C(domain), result)

def lazutkin_weight(domain, t):
    """ Returns the Lazutkin weight of domain at the point x.
    """

    return power(fprod([2, C(domain),
                        power(domain.radius_of_curv(t), fdiv(1, 3))]), -1)
