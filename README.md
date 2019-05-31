Code simulating billiard dynamics of arbitrary convex symmetric domains.

Python dependencies:
 - matplotlib
 - mpmath
 - SciPy

Visualize the billiard dynamics of an arbitrary symmetric convex domain by first specifying its Fourier coefficients in a flat text file. The nth line should correspond to a floating point representation of the coefficient of the nth cosine term. Run plot.py in Python 3 and specify the path of the text file created previously along with period desired when prompted.

To compute $`||T_{*_R} - Id||_\gamma`$, run norm_test.py. The prompt for mp.dps specifies the decimal place precision of all computations, q and j are as specified in the paper according to the operator norm (page 295), and gamma is $`\gamma`$ in the norm computation.

T_star_norm_test.py computes $`T_{*_R}`$ on a function with simple fourier coefficients (see code for details) and allows the user to compute $`|\cdot|_\gamma`$ of the result with a specified $`\gamma`$.

The basic structure of the code begins with domain.py in which the domain is fully specified, notice methods for importing Fourier coefficients, computing 2D coordinates, and various relevant derivatives.

billiards.py handles all things relevant to the billiard dynamics of a Domain, including functions to compute a singular or iterated bounce, generate orbits, and compute orbit paths.

isospectral_operators.py contains all the Isospectral operators specified in the paper as well as those modified by the Lazutkin weight.

lazutkin_coordinates.py handles all computations relevant to the Lazutkin coordinates.