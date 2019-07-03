from mpmath import *
from domain import *
from isospectral_operators import *

def fourier_basis(j):
    """ Returns the function e_j(x) = cos (2 pi j x) corresponding to the
        fourier basis.
    """

    return lambda x: cos(fprod([2, pi, j, x]))

def T_lazutkin_matrix(domain, i, j):
    """ Returns the ith element of T_lazutkin(e_j).
    """

    return l_tilde_q(domain, i, fourier_basis(j))

def general_operator_gamma_norm(matrix, gamma, max_j, max_q):
    """ Returns the gamma operator norm of matrix, summing up to max_j and
        considering the sup up to max_q. Assumed that matrix is a function
        accepting two arguments i,j and not an array () for efficiency.
    """

    max_j_sum = -1
    q = 1

    while(q < max_q):
        temp_j_sum = nsum(lambda j: fprod([power(q, gamma), power(j, -gamma),
                                           fabs(matrix(q, j))]), [1, max_j])
        max_j_sum = temp_j_sum if temp_j_sum > max_j_sum else max_j_sum
        q += 1
    return max_j_sum

def general_operator_gamma_to_sup_norm(matrix, gamma, max_j, max_q):
    """ Returns the gamma operator norm of matrix, summing up to max_j and
        considering the sup up to max_q. Assumed that matrix is a function
        accepting two arguments i,j and not an array () for efficiency.
    """

    max_j_sum = -1
    q = 1

    while(q < max_q):
        temp_j_sum = nsum(lambda j: fprod([power(j, -gamma),
                                           fabs(matrix(q, j))]), [1, max_j])
        max_j_sum = temp_j_sum if temp_j_sum > max_j_sum else max_j_sum
        q += 1
    return max_j_sum

def operator_diff_finite_sum(T_star_R_matrix, gamma, q, max_j):
    """ Returns the sum of the qth row of the matrix of T_star_R up to the
        max_jth element, scaled according to the h_gamma norm.
    """
    id = lambda i, j: 1 if i == j else 0
    sum_list = []
    j = 1

    for j in range(1, max_j):
        sum_list.append(fprod([power(q, gamma), power(j, -gamma),
                               fabs(T_star_R_matrix[j][q] - id(q, j))]))
    return sum_list

def T_star_R_operator_diff_from_id(domain, gamma, max_j, max_q):
    """ Returns the gamma operator norm T_star_R - Id, summing up to max_j and
        considering the sup up to max_q.
    """

    # Keep a dummy element for index readability, indexed by column, row
    T_star_R_matrix = [0]
    j = 1
    while (j < max_j):
        T_star_R_matrix.append(T_star_R(domain, fourier_basis(j), max_q))
        j += 1

    max_j_sum = -1
    q = 1
    while (q < max_q):

        temp_j_sum = fsum(operator_diff_finite_sum(T_star_R_matrix, gamma, q, max_j))
        max_j_sum = temp_j_sum if temp_j_sum > max_j_sum else max_j_sum
        q += 1
    return max_j_sum

def finite_h_star_gamma_norm(vector, gamma):
    """ Returns the h_star_gamma norm of vector.
    """

    norm_max = 0

    for i in range(len(vector)):
        if (fmul(power(i, gamma), abs(vector[i])) > norm_max):

            norm_max = fmul(power(i, gamma), abs(vector[i]))
    return norm_max
