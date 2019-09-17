Code simulating billiard dynamics of arbitrary convex symmetric domains.

Python dependencies:
 - matplotlib
 - mpmath
 - SciPy

Visualize the billiard dynamics of an arbitrary symmetric convex
domain by first specifying its Fourier coefficients in a flat text
file. The nth line should correspond to a floating point
representation of the coefficient of the nth cosine term. Run plot.py
in Python 3 and specify the path of the text file created previously
along with period desired when prompted.

To compute <img src="http://www.sciweavers.org/tex2img.php?eq=%7C%7C%5Cmathcal%7BT%7D_%7B%2A_R%7D%20-%20Id%7C%7C_%5Cgamma&bc=White&fc=Black&im=jpg&fs=12&ff=ccfonts,eulervm&edit=0" align="center" border="0" alt="||\mathcal{T}_{*_R} - Id||_\gamma" width="87" height="19" />, run norm_test.py. The prompt
for mp.dps specifies the decimal place precision of all computations,
q and j are as specified in the paper according to the operator norm
(page 295), and gamma is <img src="http://www.sciweavers.org/tex2img.php?eq=%5Cgamma&bc=White&fc=Black&im=jpg&fs=12&ff=ccfonts,eulervm&edit=0" align="center" border="0" alt="\gamma" width="12" height="14" /> in the norm computation.

T_star_norm_test.py computes <img src="http://www.sciweavers.org/tex2img.php?eq=%5Cmathcal%7BT%7D_%7B%2A_R%7D&bc=White&fc=Black&im=jpg&fs=12&ff=ccfonts,eulervm&edit=0" align="center" border="0" alt="\mathcal{T}_{*_R}" width="25" height="19" /> on a function with simple
fourier coefficients (see code for details) and allows the user to
compute <img src="http://www.sciweavers.org/tex2img.php?eq=%7C%5Ccdot%7C_%5Cgamma&bc=White&fc=Black&im=jpg&fs=12&ff=ccfonts,eulervm&edit=0" align="center" border="0" alt="|\cdot|_\gamma" width="32" height="19" /> of the result with a specified <img src="http://www.sciweavers.org/tex2img.php?eq=%5Cgamma&bc=White&fc=Black&im=jpg&fs=12&ff=ccfonts,eulervm&edit=0" align="center" border="0" alt="\gamma" width="12" height="14" />.

The basic structure of the code begins with domain.py in which the
domain is fully specified, notice methods for importing Fourier
coefficients, computing 2D coordinates, and various relevant
derivatives.

billiards.py handles all things relevant to the billiard dynamics of a
Domain, including functions to compute a singular or iterated bounce,
generate orbits, and compute orbit paths.

isospectral_operators.py contains all the Isospectral operators
specified in the paper as well as those modified by the Lazutkin
weight.

lazutkin_coordinates.py handles all computations relevant to the
Lazutkin coordinates.

operator_norms.py contains all functions relevant to norms of the
isospectral (and general) linear operators. Any modification of norm
computations can be made here and here alone.
