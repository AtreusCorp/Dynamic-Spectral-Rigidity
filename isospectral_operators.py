from mpmath import *
from billiards import *

def l_q(domain, q, function):
	""" Returns the result of evaluating the functional l_q on function,
	    as defined in the Dynamic Spectral Rigidity paper (TODO: SAY THIS BETTER).
	"""

	orbit = compute_q_bounce_path(domain, q, 0)
	summands = [fmul(function(point[0]), sin(point[1])) for point in orbit]
	return fsum(summands)