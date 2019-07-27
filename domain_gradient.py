from operator_norms import *

def compute_domain_grad(fourier_coeff, grad_approx_step=0.5, gamma=3.99, 
                        max_j=15, max_q=10):

    """ Return the gradient of T_star_R_operator_diff_from_id with respect
        to fourier coefficients of the domain at point fourier_coeff.
    """

    base_point = Domain()
    test_point = Domain()
    base_point.import_fourier(fourier_coeff)
    base_val = T_star_R_operator_diff_from_id(base_point, gamma, max_j, max_q)
    gradient = [0]
    i = 1

    while (i < len(fourier_coeff)):

        new_fourier = fourier_coeff.copy()
        new_fourier[i] += grad_approx_step / (10 ** i)
        test_point.import_fourier(new_fourier)
        test_val_1 = T_star_R_operator_diff_from_id(test_point, gamma, max_j,
                                                    max_q)
        new_fourier[i] -= 2 * grad_approx_step / (10 ** i)
        test_point.import_fourier(new_fourier)
        test_val_2 = T_star_R_operator_diff_from_id(test_point, gamma, max_j,
                                                    max_q)
        gradient.append((test_val_1 - test_val_2) / (2 * grad_approx_step / (10 ** i)))
        i += 1

    return gradient