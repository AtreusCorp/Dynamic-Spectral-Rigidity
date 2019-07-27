from domain import *

def c_k_domain_diff(domain_1, domain_2, k):
    """ Compute the C^k norm of domain_1.radius - domain_2.radius.
    """

    domain_1_kth_deriv = lambda s: matrix((diff(lambda t: domain_1.polar(t)[0], s, n=k), 
                                           diff(lambda t: domain_1.polar(t)[1], s, n=k)))
    domain_2_kth_deriv = lambda s: matrix((diff(lambda t: domain_2.polar(t)[0], s, n=k), 
                                           diff(lambda t: domain_2.polar(t)[1], s, n=k)))
    norm_pre_sup = lambda s: norm(domain_1_kth_deriv(s) - domain_2_kth_deriv(s))
    curvature_maximum = fdiv(1, minimize_scalar(lambda s: fdiv(1, norm_pre_sup(s)), bracket=(0, 1 / 2)).fun)
    return curvature_maximum