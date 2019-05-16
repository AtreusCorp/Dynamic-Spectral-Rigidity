Code simulating billiard dynamics of arbitrary convex symmetric domains.

Python dependencies:
 - matplotlib
 - mpmath

Visualize the billiard dynamics of an arbitrary symmetric convex domain by first specifying its Fourier coefficients in a flat text file. The nth line should correspond to a floating point representation of the coefficient of the nth cosine term. Run plot.py in Python 3 and specify the path of the text file created previously along with period desired when prompted.

The basic structure of the code begins with domain.py in which the domain is fully specified, notice methods for importing Fourier coefficients, computing 2D coordinates, and various relevant derivatives.

billiards.py handles all things relevant to the billiard dynamics of a Domain, including functions to compute a singular or iterated bounce, generate orbits, and compute orbit paths.

isospectral_operators.py contains all the Isospectral operators specified in the paper as well as those modified by the Lazutkin weight.

lazutkin_coordinates.py handles all computations relevant to the Lazutkin coordinates.