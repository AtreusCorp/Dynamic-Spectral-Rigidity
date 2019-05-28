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

    integrand_summand_1 = lambda x: fmul(fdiv(-1, 4), power(domain.radius(x), 
                                                            fdiv(-2, 3)))

    # Computing the change of coordinates for the first derivatve
    radius_derivative_lazut = lambda x: fdiv(domain.radius_derivative(x), 
        lazutkin_param_non_arc_deriv(domain, lazutkin_param_non_arc(domain, x)))

    # Computing the change of coordinates for the second derivative
    radius_second_deriv_lazut_term_1 = lambda x: fdiv(domain.radius_second_derivative(x), 
        power(lazutkin_param_non_arc_deriv(domain, lazutkin_param_non_arc(domain, x)), 2))

    radius_second_deriv_lazut_term_2 = lambda x: fmul(fdiv(domain.radius_derivative(x), 
        power(lazutkin_param_non_arc_deriv(domain, lazutkin_param_non_arc(domain, x)), 2)),
    lazutkin_param_non_arc_second_deriv(domain, lazutkin_param_non_arc(domain, x)))

    radius_second_derivative_lazutkin = lambda x: fsub(radius_second_deriv_lazut_term_1(x),
                                                       radius_second_deriv_lazut_term_2(x))
    integrand_summand_2 = lambda x: fprod([fdiv(1, 6), 
                                           power(domain.radius(x), fdiv(1, 3)), 
                                           radius_second_derivative_lazutkin(x)])

    integrand_summand_3 = lambda x: fprod([fdiv(-1, 9), 
        power(domain.radius(x), fdiv(-2, 3)), 
        power(radius_derivative_lazut(x), 2)])
    integrand_sum = lambda x: fsum([integrand_summand_1(x), 
                                    integrand_summand_2(x),
                                    integrand_summand_3(x)])
    integrand = lambda x: fprod([integrand_sum(x), 
                                 function(lazutkin_param_non_arc(domain, x)), 
                                 lazutkin_param_non_arc_deriv(domain, x)])
    return fmul(fdiv(1, fmul(6, power(C(domain), 2))), quad(integrand, [power(10, -2 * mp.dps), 1]))

def T_star_R(domain, function, precision):
    """ Returns a list of length precision given by a component of the 
        linearized isospectral operator modified by the Lazutkin weight.
    """

    T_val = T_lazutkin(domain, function, precision)
    l_0 = l_q_lazutkin(domain, 0, function)
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