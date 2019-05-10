from mpmath import *
from billiards import *

def l_q(domain, q, function):
	""" Returns the result of evaluating the functional l_q on function,
	    as defined in the Dynamic Spectral Rigidity paper (TODO: SAY THIS BETTER).
	"""
	

	if (q == 0):
		integrand = lambda x: fmul(fdiv(function(x), domain.radius(x)), 
								   norm(domain.polar_gradient(x)))

		# Here we make the assumption that the length of the boundary of
		# domain is 1
		return quad(integrand, [0, 1])

	elif (q == 1):
		return function(0)

	orbit = compute_q_bounce_path(domain, q, 0)
	summands = [fmul(function(point[0]), sin(point[1])) for point in orbit]
	return fsum(summands)

def T(domain, function, precision):
	""" Returns a list of length precision given by the linearized isospectral 
		operator.
	"""

	output = []
	q = 0

	while (q < precision):
		output.append(l_q(domain, q, function))
		
	return output
