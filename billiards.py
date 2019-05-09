from mpmath import *
from domain import Domain
from math import isinf

def bounce(domain, inc_theta, inc_angle):
    """ Returns the next point after an edge collision at inc_angle, the 
        angle of incidence inc_angle in domain. This is essentially the 
        billiard map without the cos of the angle. 
    """

    inc_angle = fmod(inc_angle, pi)

    incident_point = domain.polar(inc_theta)
    path_vector_x = fadd(fmul(cos(inc_angle), 
                         domain.polar_gradient(inc_theta)[0]),
                         fneg(fmul(sin(inc_angle), 
                         domain.polar_gradient(inc_theta)[1])))
    path_vector_y = fadd(fmul(sin(inc_angle), 
                         domain.polar_gradient(inc_theta)[0]),
                         fmul(cos(inc_angle), 
                         domain.polar_gradient(inc_theta)[1]))
    path_vector = (path_vector_x, path_vector_y)

    difference_fnc = lambda theta, t: (fsub(fadd(incident_point[0], 
                                                 fmul(t, path_vector_x)),
                                            domain.polar(theta)[0]),
                                       fsub(fadd(incident_point[1], 
                                                 fmul(t, path_vector_y)),
                                            domain.polar(theta)[1]))
    newton_start_point = [pi, 1]
    next_point = [0, 0]

    # Make sure the root found isn't the point of incidence.
    # Use density of irrational rotations of the circle to pick another
    # starting point
    while almosteq(next_point[1], 0):
        newton_start_point[0] += 1
        try:
            next_point = findroot(difference_fnc, newton_start_point)
        except:
            next_point = [0, 0]

    next_inc_theta = fmod(next_point[0], fmul(2, pi))
    next_point_t = next_point[1]

    path_vector_end = (fadd(incident_point[0], 
                            fmul(next_point_t, path_vector_x)),
                       fadd(incident_point[1], 
                            fmul(next_point_t, path_vector_y)))
    # Cosine law
    next_point_inc_angle = acos(fdiv(fdot(path_vector, 
                                          domain.polar_gradient(next_inc_theta)),
                                     fmul(norm(path_vector), 
                                          norm(domain.polar_gradient(next_inc_theta)))))
    return (next_inc_theta, next_point_inc_angle)

def bounce_q_times(domain, q, inc_theta, inc_angle):
    """ Return the value obtained by iteratively applying bounce q times. 
    """

    cur_pt = (inc_theta, inc_angle)
    i = 0
    while i < q:
        cur_pt = bounce(domain, cur_pt[0], cur_pt[1])
        i += 1
    return cur_pt

def bounce_q_times_top_half(domain, q, inc_theta, inc_angle):
    """ Return the value obtained by iteratively applying bounce q times. 
    """

    cur_pt = (inc_theta, inc_angle)
    penalty = 0
    i = 0
    while i < q:
        cur_pt = bounce(domain, cur_pt[0], cur_pt[1])

        # Incompatible bounce
        if (pi < fabs(fsub(inc_theta, cur_pt[0]))):
            penalty = fadd(penalty, fsub(fabs(fsub(inc_theta, cur_pt[0])), pi))

        i += 1
    return (cur_pt, penalty)

def check_q_bounce_angle(domain, q, inc_theta, inc_angle):
    """ Returns the difference between inc_theta after the qth and q + 1 st 
        bounces along with the penalty incurred. This is a helper to generate_orbit_odd.
    """

    last_bounce = bounce_q_times_top_half(domain, q, inc_theta, inc_angle)    
    following_bounce = bounce(domain, last_bounce[0][0], last_bounce[0][1])

    return (fsub(fadd(last_bounce[0][0], following_bounce[0]), fmul(2, pi)), last_bounce[1])



def generate_orbit_odd(domain, q, start_theta):
    """ Returns an angle theta such that bounce is periodic of period q
        on input theta, beginning at the marked point. Assumed that q is 
        odd.
    """

    bounce_end_theta_check = lambda bounce_angle: check_q_bounce_angle(domain, (q - 1) / 2, start_theta, bounce_angle)
    newton_start_point = 0
    found_root = [0]

    # Make sure the root found isn't the point of incidence.
    # Use density of irrational rotations of the circle to pick another
    # starting point
    while almosteq(found_root[0], start_theta) or abs(fsub(found_root[0], start_theta)) >= fdiv(pi, 2):
        newton_start_point += 1 / 2
        try:
            found_root = findroot(bounce_end_theta_check, newton_start_point)
        except:
            found_root[0] = 0
    return found_root[0]

def check_final_bounce_angle(domain, q, inc_theta, inc_angle):
    """ Returns the penalty and final polar angle after q bounces starting from 
        inc_theta and inc_angle with 2pi subtracted. Helper function to 
        generate_orbit_even.
    """

    final_theta = bounce_q_times_top_half(domain, q, inc_theta, inc_angle)

    return (fsub(final_theta[0][0], fadd(inc_theta, pi)), final_theta[1])

def generate_orbit_even(domain, q, start_theta):
    """ Returns an angle theta such that bounce is periodic of period q
        on input theta, beginning at the marked point. Assumed that q is 
        even.
    """

    if q == 2:
        return fdiv(pi, 2)

    bounce_end = lambda bounce_angle: check_final_bounce_angle(domain, 
        q / 2, start_theta, bounce_angle)
    newton_start_point = 0
    found_root = [0]

    # Make sure the root found isn't the point of incidence.
    # Use density of irrational rotations of the circle to pick another
    # starting point
    while almosteq(found_root[0], start_theta) or abs(fsub(found_root[0], start_theta)) >= fdiv(pi, 2):
        newton_start_point += 1 / 2
        try:
            found_root = findroot(bounce_end, newton_start_point)
        except:
            found_root[0] = 0
    return found_root[0]

def generate_orbit(domain, q, start_theta):
    """ Returns an angle theta such that bounce is periodic of period q
        on input theta, beginning at the marked point.
    """

    if (q % 2 == 0):
        return generate_orbit_even(domain, q, start_theta)
    else:
        return generate_orbit_odd(domain, q, start_theta)