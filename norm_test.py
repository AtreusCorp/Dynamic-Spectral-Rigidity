from operator_norms import *

mp.dps = int(input("Enter mp.dps: "))
domain = Domain()
fourier_path = input("Enter a path to a Fourier coefficient file: ")
domain.import_fourier_from_file(fourier_path)
max_j = int(input("Enter max j: "))
max_q = int(input("Enter max q: "))
gamma = float(input("Enter gamma: "))
print(T_star_R_operator_diff_from_id(domain, gamma, max_j, max_q))