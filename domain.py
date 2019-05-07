from mpmath import *

#TODO specify precision
mp.prec = 53                #[default: 53]
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

        #TODO Precision
        series_terms = [fmul(self.fourier[n], cos(fmul(n, theta))) 
                        for n in range(0, len(self.fourier))]
        return fsum(series_terms)

    def radius_derivative(self, theta):
        """ Computes the derivative of the radius function at a given angle. 
            Assumed that import_fourier has been called. 
        """

        series_terms = [fprod([- n, self.fourier[n], sin(fmul(n, theta))]) 
                        for n in range(0, len(self.fourier))]
        return fsum(series_terms)

    def polar(self, theta):
        """ Provides an interface for polar coordinates. 
            Assumed that import_fourier has been called.
        """

        return (fmul(self.radius(theta), cos(theta)),  
                fmul(self.radius(theta), sin(theta)))

    def polar_gradient(self, theta):
        """ Returns the gradient of the polar map with respect to theta.
        """

        r = self.radius(theta)
        r_prime = self.radius_derivative(theta)
        grad_x = fsub(fmul(r_prime, cos(theta)), fmul(r, sin(theta)))
        grad_y = fadd(fmul(r_prime, sin(theta)), fmul(r, cos(theta)))

        return (grad_x, grad_y)