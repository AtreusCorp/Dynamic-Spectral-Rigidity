from mpmath import *

#TODO specify precision
mp.dps = 15                 #[default: 15]

class Domain:
    """ Assumed to be a Z_2 symmetric convex domain.
    """
    def __init__(self):

        self.fourier = []

    def import_fourier(self, path):
        """ Parses a textfile located at path containing fourier coefficients 
            for the domain. Coefficients are expected one per line. 
        """

        file = open(path, "r")

        for line in file:
            self.fourier.append(mpf(line))
        return

    def radius(self, theta):
        """ Computes the radius function at a given angle. Assumed that
            import_fourier has been called. 
        """

        return fourierval((self.fourier, [0]), [0, 1], theta)

    def radius_derivative(self, theta):
        """ Computes the derivative of the radius function at a given angle. 
            Assumed that import_fourier has been called. 
        """

        series_terms = [fprod([-n, fmul(2, pi), self.fourier[n]]) 
                        for n in range(0, len(self.fourier))]
        return fourierval(([0], series_terms), [0, 1], theta)

    def polar(self, theta):
        """ Provides an interface for polar coordinates. 
            Assumed that import_fourier has been called.
        """
        theta_rescaled = fprod([theta, 2, pi])
        return (fmul(self.radius(theta), cos(theta_rescaled)),  
                fmul(self.radius(theta), sin(theta_rescaled)))

    def polar_gradient(self, theta):
        """ Returns the gradient of the polar map with respect to theta.
        """

        r = self.radius(theta)
        r_prime = self.radius_derivative(theta)
        theta_rescaled = fprod([theta, 2, pi])
        grad_x = fsub(fmul(r_prime, cos(theta_rescaled)), fprod([r, 2, pi, sin(theta_rescaled)]))
        grad_y = fadd(fmul(r_prime, sin(theta_rescaled)), fprod([r, 2, pi, cos(theta_rescaled)]))

        return (grad_x, grad_y)

def arc_length_coords(domain, x):
    """ Returns the transformation taking x (in [0, 2 pi]) to the 
        corresponding point in the parameterization by arc length.
    """

    return quad(lambda t: norm(domain.polar_gradient(t)), [0, x])