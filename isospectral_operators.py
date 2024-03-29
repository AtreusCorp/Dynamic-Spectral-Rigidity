from mpmath import *
from domain import *
from billiards import *
from lazutkin_coordinates import *

def l_q(domain, q, function):
    """ Returns the result of evaluating the functional l_q. (TODO: SAY THIS BETTER).
    """

    if (q == 0):
        integrand = lambda x: fprod([fdiv(function(x), domain.radius_of_curv(x)),
                                     norm(domain.polar_gradient(x))])

        return quad(integrand, [0, 1])

    elif (q == 1):
        return function(0)

    orbit = generate_orbit(domain, q)
    summands = [fmul(function(point[0]), sin(point[1])) for point in orbit]
    return fsum(summands)

def l_tilde_q(domain, q, function):
    """ Returns the output of the linearized map modified by the Lazutkin
        weight.
    """
    if (q == 0):
        integrand = lambda t: fmul(2, function(t))
        return quad(integrand, [0, 1])

    fun_modified = lambda t: fdiv(function(lazutkin_param_non_arc(domain, t)),
                                  lazutkin_weight(domain, t))

    return l_q(domain, q, fun_modified)

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
        output.append(l_tilde_q(domain, q, function))
        q += 1
    return output

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

def l_bullet(domain, function):
    """ Returns the value of ell_bullet on function. TODO: Make sure the formula
        is correct.
    """

    integrand_summand_1 = lambda t: fmul(fdiv(-1, fmul(24, power(C(domain), 2))), 
                                         power(domain.radius_of_curv(t), fdiv(-2, 3)))

    # Computing the change of coordinates for the first derivatve
    radius_derivative_lazut = lambda t: fdiv(domain.radius_of_curv_deriv(t), 
                                             lazutkin_param_non_arc_deriv(domain, t))

    # Computing the change of coordinates for the second derivative
    radius_second_derivative_lazutkin_1 = lambda t: fsub(domain.radius_of_curv_second_deriv(t), 
                                                         fmul(radius_derivative_lazut(t), 
                                                              lazutkin_param_non_arc_second_deriv(domain, t)))
    radius_second_derivative_lazutkin = lambda t: fdiv(radius_second_derivative_lazutkin_1(t), 
        power(lazutkin_param_non_arc_deriv(domain, t), 2))

    integrand_summand_2 = lambda t: fprod([fdiv(1, 36),
                                           power(domain.radius_of_curv(t), -1),
                                           radius_second_derivative_lazutkin(t)])

    integrand_summand_3 = lambda t: fprod([fdiv(-1, 27),
        power(domain.radius_of_curv(t), -2),
        power(radius_derivative_lazut(t), 2)])
    integrand_sum = lambda t: fsum([integrand_summand_1(t),
                                    integrand_summand_2(t),
                                    integrand_summand_3(t)])
    integrand = lambda t: fprod([integrand_sum(t),
                                 function(lazutkin_param_non_arc(domain, t)),
                                 lazutkin_param_non_arc_deriv(domain, t)])
    return  quadgl(integrand, [0, 1])

def T_star_R(domain, function, precision):
    """ Returns a list of length precision given by a component of the
        linearized isospectral operator modified by the Lazutkin weight.
    """

    T_val = T_lazutkin(domain, function, precision)
    l_0 = T_val[0]
    b_0 = T_lazutkin(domain, lambda x: 1, precision)

    # Attention: There is a small typo in the original paper in Lemma 5.3
    # b_l l_0 should have a 1/2 coefficient which is reflected here

    l_0_term = [fprod([0.5, b_term, l_0]) for b_term in b_0]
    l_bullet_val = l_bullet(domain, P_star(function))
    b_bullet_val = b_bullet(precision)
    l_bullet_term = [fmul(b_term, l_bullet_val) for b_term in b_bullet_val]
    sub_terms = [fadd(zero_term, bullet_term) for zero_term, bullet_term
        in zip(l_0_term, l_bullet_term)]
    return [fsub(t_term, sub_term) for t_term, sub_term in zip(T_val, sub_terms)]