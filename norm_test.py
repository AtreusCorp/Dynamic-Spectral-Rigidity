from isospectral_operators import *

mp.dps = int(input("Enter mp.dps: "))
domain = Domain()
fourer_path = input("Enter a path to a Fourier coefficient file: ")
domain.import_fourier(fourer_path)
max_j = int(input("Enter max j: "))
max_q = int(input("Enter max q: "))
gamma = float(input("Enter gamma: "))
print(T_star_R_operator_diff_from_id(domain, gamma, max_j, max_q))