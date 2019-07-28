from operator_norms import *

def compute_domain_grad(fourier_coeff, grad_approx_step=0.005, gamma=3.99, 
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

def euler_domain_grad(fourier_coeff, grad_approx_step=0.005, euler_step=0.0001,
                      gamma=3.99, max_j=15, max_q=10, num_steps=10):
    """ Perform Eulers method on the fourier coefficients of the domain at 
        fourier_coeff using compute_domain_grad, print each step.
    """

    base_point = Domain()
    test_point = Domain()
    test_coeff = fourier_coeff
    base_point.import_fourier(fourier_coeff)
    base_val = T_star_R_operator_diff_from_id(base_point, gamma, max_j, max_q)
    i = 0

    while (i < num_steps):

        gradient = compute_domain_grad(test_coeff, grad_approx_step, gamma, 
                                       max_j, max_q)
        for j in range(len(test_coeff)):
            test_coeff[j] -= fmul(euler_step, gradient[j])

        test_point.import_fourier(test_coeff)
        print("Point: {}, Value: {}\n".format(test_coeff, 
            T_star_R_operator_diff_from_id(test_point, gamma, max_j, max_q)))
        i += 1

def euler_domain_grad_increase(fourier_coeff, grad_approx_step=0.005, 
                               euler_step=0.0001, gamma=3.99, max_j=15, 
                               max_q=10, num_steps=10):
    """ Perform Eulers method on the fourier coefficients of the domain at 
        fourier_coeff using compute_domain_grad, print each step.
    """

    base_point = Domain()
    test_point = Domain()
    test_coeff = fourier_coeff
    base_point.import_fourier(fourier_coeff)
    new_norm = T_star_R_operator_diff_from_id(base_point, gamma, max_j, max_q)
    i = 0

    while (i < num_steps):

        gradient = compute_domain_grad(test_coeff, grad_approx_step, gamma, 
                                       max_j, max_q)
        if new_norm > 1:
            print("Method exceeded 1. Try again with new constants.\n")
            return
        elif (1 - new_norm < euler_step)
            euler_step = (1 - new_norm) * euler_step

        for j in range(len(test_coeff)):
            test_coeff[j] += fmul(euler_step, gradient[j])

        test_point.import_fourier(test_coeff)
        new_norm = T_star_R_operator_diff_from_id(test_point, gamma, max_j, 
                                                  max_q)
        print("Point: {}, Value: {}\n".format(test_coeff, new_norm))
        i += 1