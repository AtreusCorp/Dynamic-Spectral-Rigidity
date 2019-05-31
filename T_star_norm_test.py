from isospectral_operators import *
from operator_norms import *

mp.dps = int(input("Enter mp.dps: "))
domain = Domain()
fourer_path = input("Enter a path to a Fourier coefficient file: ")
domain.import_fourier(fourer_path)
function_len = int(input("Enter length of test function: "))
list = []
for i in range(1, function_len):
    list.append(pi / i**4)

function = lambda x: fourierval((list, [0]), [0, 1], x)
t_star_prec = int(input("Enter T_star_R precision: "))
data = T_star_R(domain, function, t_star_prec)
print(data)
gamma = float(input("Enter gamma: "))
print(finite_h_star_gamma_norm(data, gamma))