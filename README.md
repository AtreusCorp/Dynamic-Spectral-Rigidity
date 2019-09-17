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

To compute `||T_{*_R} - Id||_\gamma`, run norm_test.py. The prompt
for mp.dps specifies the decimal place precision of all computations,
q and j are as specified in the paper according to the operator norm
(page 295), and gamma is `\gamma` in the norm computation.

T_star_norm_test.py computes `T_{*_R}` on a function with simple
fourier coefficients (see code for details) and allows the user to
compute `|\cdot|_\gamma` of the result with a specified `\gamma`.

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

Copyright 2019 University of Toronto

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
