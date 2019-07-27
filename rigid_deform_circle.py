from operator_norms import *

mp.dps = 20
precision_step = 1
fourier_length = 20
fourier_list = [1] + [0] * fourier_length
max_j = 10
max_q = 10
gamma = 3.999999
lower_domain = Domain()
upper_domain = Domain()

file = open('test_file', 'w+')
file.write('[1, 1]\n')
i = 1
while (i < fourier_length):

    cur_min = precision_step / (10 ** i)
    cur_max = -precision_step / (10 ** i)
    upper_list = fourier_list.copy()
    lower_list = fourier_list.copy()
    lower_domain.import_fourier(lower_list)
    upper_domain.import_fourier(upper_list)
    done = 0

    while (not done):

        lower_norm = T_star_R_operator_diff_from_id(lower_domain, gamma, max_j, max_q)
        upper_norm = T_star_R_operator_diff_from_id(upper_domain, gamma, max_j, max_q)

        if (lower_norm < 1):
            cur_min -= precision_step / (10 ** i)
            lower_list[i] = cur_min
            lower_domain.import_fourier(lower_list)

        if (upper_norm < 1):
            cur_max += precision_step / (10 ** i)
            upper_list[i] = cur_max
            upper_domain.import_fourier(upper_list)

        if (upper_norm >= 1 and lower_norm >= 1):
            done = 1

    file.write('({}, {})\n'.format(cur_min, cur_max))
    file.flush()
    i += 1