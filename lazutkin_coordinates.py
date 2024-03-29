from mpmath import *
from domain import *

def C(domain):
    """ Returns the C coefficient for the Lazutkin paramterization.
    """
    if domain.lazutkin_C:
        return domain.lazutkin_C

    integrand = lambda t: fmul(power(domain.radius_of_curv(t), fdiv(-2, 3)),
                               norm(domain.polar_gradient(t)))
    domain.lazutkin_C = power(quadgl(integrand, [0, 1]), -1)
    return domain.lazutkin_C

def lazutkin_param_non_arc(domain, t):
    """ Returns the Lazutkin parameter for domain corresponding to the
        polar angle t.

    :param domain: the domain
    :param t: polar angle of the point under consideration

    """

    if (t == 0 or t == 1):
        return t

    closest_mesh_edge = min(domain.lazutkin_cache.keys(), 
                            key=lambda elt: abs(elt - t))
    integral = 0
    integral += domain.lazutkin_cache[closest_mesh_edge]

    if (almosteq(closest_mesh_edge, t)):
        return integral

    integrand = lambda s: fmul(power(domain.radius_of_curv(s), fdiv(-2, 3)),
                               norm(domain.polar_gradient(s)))
    return fadd(integral, 
                fmul(C(domain), quadgl(integrand, [closest_mesh_edge, t])))

def inv_lazutkin_param_non_arc(domain, x):
    """ Returns the inverse of the transformation taking x (in [0, 2 pi]) to the 
        corresponding point in the parameterization by arc length.
    """

    closest_mesh_edge = min(domain.lazutkin_cache.keys(), 
                            key=lambda elt: abs(domain.lazutkin_cache[elt] - x))
    if (x < domain.lazutkin_cache[closest_mesh_edge]):
        bounds = [closest_mesh_edge - domain.lazutkin_mesh, closest_mesh_edge]
    elif (x > domain.lazutkin_cache[closest_mesh_edge]):
        bounds = [closest_mesh_edge, closest_mesh_edge + domain.lazutkin_mesh]

    # x == domain.lazutkin_cache[closest_mesh_edge]
    else:
        return closest_mesh_edge

    return findroot(lambda t: lazutkin_param_non_arc(domain, t) - x, bounds, solver='illinois')

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