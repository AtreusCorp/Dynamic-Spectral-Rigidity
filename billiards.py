from mpmath import *
from domain import Domain
from lazutkin_coordinates import *
from scipy.optimize import minimize_scalar

def bounce_compare(domain, path_vector, incident_point, t_0):
    """ A helper for bounce which checks the validity of a bounce point
        candidate. Return True if halving t_0 brings one closer to
        the domain.
    """

    polar_1 = lambda theta: domain.polar(theta)[0]
    polar_2 = lambda theta: domain.polar(theta)[1]
    path_vec_1 = lambda t: fadd(fmul(path_vector[0], t),
                                incident_point[0])
    path_vec_2 = lambda t: fadd(fmul(path_vector[1], t),
                                incident_point[1])
    G = lambda theta: fadd(power(fsub(polar_1(theta), path_vec_1(t_0)), 2),
                           power(fsub(polar_2(theta), path_vec_2(t_0)), 2))
    F = lambda theta: fadd(power(fsub(polar_1(theta), path_vec_1(fdiv(t_0, 2))), 2),
                           power(fsub(polar_2(theta), path_vec_2(fdiv(t_0, 2))), 2))
    return minimize_scalar(G).fun > minimize_scalar(F).fun

def bounce(domain, inc_theta, inc_angle, guess_theta = 1 / 2):
    """ Returns the next point after an edge collision at inc_angle, the
        angle of incidence inc_angle in domain. This is essentially the
        billiard map without the cos of the angle.
    """

    inc_angle = fmod(inc_angle, pi)
    if (almosteq(inc_angle, 0)):
        return (inc_theta, inc_angle)

    incident_point = domain.polar(inc_theta)
    incident_grad = domain.polar_gradient(inc_theta)
    path_vector_x = fadd(fmul(cos(inc_angle), incident_grad[0]),
                         fneg(fmul(sin(inc_angle), incident_grad[1])))
    path_vector_y = fadd(fmul(sin(inc_angle), incident_grad[0]),
                         fmul(cos(inc_angle), incident_grad[1]))
    path_vector = (path_vector_x, path_vector_y)
    difference_fnc = lambda theta, t: (fsub(fadd(incident_point[0],
                                                 fmul(t, path_vector_x)),
                                            domain.polar(theta)[0]),
                                       fsub(fadd(incident_point[1],
                                                 fmul(t, path_vector_y)),
                                            domain.polar(theta)[1]))
    newton_start_point = [guess_theta, 1]
    next_point = [0, 0]

    # Make sure the root found isn't the point of incidence.
    # Use density of irrational rotations of the circle to pick another
    # starting point
    while (almosteq(next_point[1], 0)
            or almosteq(fmod(next_point[0], 1), inc_theta)
            or bounce_compare(domain, path_vector, incident_point, next_point[1])):
        try:
            next_point = findroot(difference_fnc, newton_start_point,
                                  tol = power(10, -2 * mp.dps))
            next_point[0] = fmod(next_point[0], 1)
        except:
            next_point = [0, 0]
        finally:
            newton_start_point[0] += fdiv(pi, 8)

    next_inc_theta = fmod(next_point[0], 1)
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

def bounce_q_times(domain, q, inc_angle):
    """ Return the value obtained by iteratively applying bounce q times,
        beginning from angle 0.
    """

    cur_pt = (0, inc_angle)
    i = 0
    while i < q:
        cur_pt = bounce(domain, cur_pt[0], cur_pt[1])
        i += 1
    return cur_pt

def bounce_q_times_top_half(domain, q, inc_angle):
    """ Return the value obtained by iteratively applying bounce q times,
        starting from angle 0 and adding a penalty for bounces that land in
        the lower half of the domain.
    """

    cur_pt = (0, inc_angle)
    penalty = 0
    i = 0
    while i < q:
        cur_pt = bounce(domain, cur_pt[0], cur_pt[1])

        # Incompatible bounce
        if ((1 / 2) < fabs(cur_pt[0])):
            penalty = fadd(penalty, fsub(fabs(cur_pt[0]), 1 / 2))

        i += 1
    return (cur_pt, penalty)

def check_q_bounce_down(domain, q, inc_angle):
    """ Returns the difference in polar angle between the q and q + 1 st bounces
        along with the penalty incurred by bouncing below the top half of the
        domain.
    """

    last_bounce = bounce_q_times_top_half(domain, q, inc_angle)
    following_bounce = bounce(domain, last_bounce[0][0], last_bounce[0][1])
    return (fsub(fadd(last_bounce[0][0], following_bounce[0]), 1),
                 last_bounce[1])

def generate_orbit_odd(domain, q):
    """ Returns an angle theta such that bounce is periodic of period q
        on input theta, beginning at the marked point. Assumed that q is
        odd.
    """

    bounce_end = lambda bounce_angle: check_q_bounce_down(domain, (q - 1) / 2,
                                                          bounce_angle)
    newton_start_point = lazutkin_weight(domain, 0) / q
    found_root = 0

    # Make sure the root found isn't the point of incidence.
    # Use density of irrational rotations of the circle to pick another
    # starting point
    while almosteq(found_root, 0) or abs(found_root) >= fdiv(pi, 2):
        try:
            found_root = fmod(findroot(bounce_end, newton_start_point)[0], pi)
        except:
            found_root = 0
        finally:
            newton_start_point += 1 / 2
            newton_start_point = fmod(newton_start_point, fdiv(pi, 2))
    return found_root

def check_q_bounce_opp(domain, q, inc_angle):
    """ Returns the difference between the polar angle of the qth bounce and 1 / 2
        along with the penalty incurred from bouncing below the top half of the
        domain. Helper function to generate_orbit_even.
    """

    final_theta = bounce_q_times_top_half(domain, q, inc_angle)
    return (fsub(final_theta[0][0], 1 / 2), final_theta[1])

def generate_orbit_even(domain, q):
    """ Returns an angle theta such that bounce is periodic of period q
        on input theta, beginning at the marked point. Assumed that q is
        even.
    """

    if q == 2:
        return fdiv(pi, 2)

    bounce_end = lambda bounce_angle: check_q_bounce_opp(domain, q / 2,
                                                         bounce_angle)
    newton_start_point = lazutkin_weight(domain, 0) / q
    found_root = 0

    # Make sure the root found isn't the point of incidence.
    # Use density of irrational rotations of the circle to pick another
    # starting point
    while almosteq(found_root, 0) or abs(found_root) >= fdiv(pi, 2):
        try:
            found_root = fmod(findroot(bounce_end, newton_start_point)[0], pi)
        except:
            found_root = 0
        finally:
            newton_start_point += 1 / 2
            newton_start_point = fmod(newton_start_point, fdiv(pi, 2))
    return found_root

def compute_q_bounce_path(domain, q):
    """ Returns a list of bounce points for an orbit with period q.
    """

    if (q % 2 == 0):
        generated_start_bounce = generate_orbit_even(domain, q)
    else:
        generated_start_bounce = generate_orbit_odd(domain, q)

    cur_bounce = [0, generated_start_bounce]
    orbit = [(cur_bounce[0], cur_bounce[1])]
    i = 0

    while(i < q - 1):
        cur_bounce = bounce(domain, cur_bounce[0], cur_bounce[1])
        orbit.append((cur_bounce[0], cur_bounce[1]))
        i += 1

    return orbit

def length_gradient(domain, i, s):
    """ Returns the ith partial derivative (with respect to t) of the length
        function for the qth orbit at the point s.
    """

    if (i >= len(s)):
        raise Exception("i must be < len(s).\n")

    if (i == 0):
        return 0

    if (i == len(s) - 1):

        # Add a dummy value to allow for odd case
        if (not almosteq(s[i], (1 / 2))):
            s = s.copy()
            s.append(1 - s[i])
        else:
            return 0

    point_grad = domain.polar_gradient(s[i])
    point_i_sub_1 = domain.polar(s[i - 1])
    point_i = domain.polar(s[i])
    point_i_plus_1 = domain.polar(s[i + 1])
    diff_vector_incoming = (point_i[0] - point_i_sub_1[0],
                            point_i[1] - point_i_sub_1[1])
    diff_vector_outgoing = (point_i_plus_1[0] - point_i[0],
                            point_i_plus_1[1] - point_i[1])
    incident_angle = fdiv(fdot(point_grad, diff_vector_incoming),
                          fmul(norm(point_grad), norm(diff_vector_incoming)))
    outgoing_angle = fdiv(fdot(point_grad, diff_vector_outgoing),
                          fmul(norm(point_grad), norm(diff_vector_outgoing)))
    return fmul(norm(domain.polar_gradient(s[i])),
                fsub(incident_angle, outgoing_angle))

def cos_outgoing_angle(domain, i, s):
    """ Returns cos of the outgoing angle at position i of the orbit s.
    """

    if (i >= len(s)):
        raise Exception("i must be < len(s).\n")

    if (i == len(s) - 1):

        # Add a dummy value to allow for odd and even cases
        # A slight hack to make the function more general
        s = s.copy()

        if (not almosteq(s[i], (1 / 2))):
            s.append(1 - s[i])
        else:
            s.append(1 - s[i - 1])

    point_grad = domain.polar_gradient(s[i])
    point_i = domain.polar(s[i])
    point_i_plus_1 = domain.polar(s[i + 1])
    diff_vector_outgoing = (point_i_plus_1[0] - point_i[0],
                            point_i_plus_1[1] - point_i[1])
    outgoing_angle = fdiv(fdot(point_grad, diff_vector_outgoing),
                          fmul(norm(point_grad), norm(diff_vector_outgoing)))
    return outgoing_angle

def compute_q_bounce_path_euler(domain, q, epsilon):
    """ Returns a list of bounce points for an orbit with period q,
        uses the euler method to optimize the length function.
    """

    cur_point = [inv_lazutkin_param_non_arc(domain, fdiv(i, q))
                 for i in range(q // 2 + 1)]
    cur_point[-1] = 1 / 2 if (q % 2 == 0) else cur_point[-1]
    grad = [length_gradient(domain, i, cur_point) for i in range(len(cur_point))]
    grad_sup_norm = norm(grad, p=inf)

    while (grad_sup_norm > power(10, -(mp.dps - 1))):

        new_point = []

        for i in range(len(cur_point)):
            new_point[i].append(fadd(cur_point[i],
                fmul(epsilon, length_gradient(domain, i, cur_point))))

        cur_point = new_point
        grad = [length_gradient(domain, i, cur_point) for i in range(len(cur_point))]
        grad_sup_norm = norm(grad, p=inf)

    # Lots of list tricks to construct the bounce path given the information
    # found by euler's method. Relies heavily on symmetry
    outgoing_angles = [acos(cos_outgoing_angle(domain, i, cur_point))
        for i in range(len(cur_point))]
    cur_point = list(zip(cur_point, outgoing_angles))

    if (q % 2):
        return_point = cur_point + [(1 - point[0], point[1]) for point in reversed(cur_point[1:])]
    else:
        return_point = cur_point + [(1 - point[0], point[1]) for point in reversed(cur_point[1:-1])]

    # Ensure that the computed path has winding number 1
    if (return_point != sorted(return_point, key=lambda x: x[0])):
            raise Exception("Euler method converged to bad orbit.")
            
    return return_point

def generate_orbit(domain, q, method='euler', epsilon=0.1, euler_try_count=3):
    """ Return an angle theta such that bounce is periodic of period q for 
        domain on input theta, beginning at the marked point.

        method -- the method used to compute bounce paths (default 'euler')
        epsilon -- if using the euler method, this is the coefficient for the
                   gradient step (default 0.1)
        euler_try_count -- the number of times to retry using the euler method
                           before defaulting to a less efficient method 
                           (default 3)
    """

    if (q in domain.orbits.keys()):
        return domain.orbits[q]

    if method == 'euler':
        done = 0
        
        # Try to compute the bounce path using euler's method euler_try_count
        # times. If no success, default to the original slower method
        while (done < 1):

            if (done <= -euler_try_count):
                orbit = compute_q_bounce_path(domain, q)
                done = 1
            else:
                try:
                    orbit = compute_q_bounce_path_euler(domain, q, epsilon)
                    done = 1
                except:
                    epsilon /= 2
                    done -= 1
        domain.orbits[q] = orbit
    else:
        domain.orbits[q] = compute_q_bounce_path(domain, q)
    return domain.orbits[q]