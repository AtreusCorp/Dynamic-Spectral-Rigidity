from mpmath import *
from lazutkin_coordinates import *
from scipy.optimize import minimize_scalar

# Specify precision
mp.dps = 15                 #[default: 15]
curvature_lower_bound = 0.1
curvature_upper_bound = 10000
lazutkin_mesh = 0.05

def arc_length_coords(domain, x):
    """ Returns the transformation taking x (in [0, 1]) to the 
        corresponding point in the parameterization by arc length.
    """

    return quad(lambda t: norm(domain.polar_gradient(t)), [0, x])

def inv_arc_length_coords(domain, s):
    """ Returns the inverse of the transformation taking x (in [0, 2 pi]) to 
        the corresponding point in the parameterization by arc length.
    """

    return findroot(lambda x: arc_length_coords(domain, x) - s, 0.5)

class Domain:
    """ Assumed to be a Z_2 symmetric convex domain.
    """
    def __init__(self):

        self.fourier = []
        self.lazutkin_cache = {}
        self.orbits = {}

    def _clean_domain(self):
        """ Clears the instance variables for self.
        """

        self.fourier = []
        self.lazutkin_cache = {}
        self.orbits = {}

    def _cache_lazutkin(self):
        """ Caches lazutkin values according to the lazutkin_mesh 
            parameter into domain. Helper to import_fourier.
        """

        self.lazutkin_cache[0] = 0
        cur_val = lazutkin_mesh

        # TODO: Make this smarter
        while (cur_val <= 1):
            self.lazutkin_cache[cur_val] = lazutkin_param_non_arc(self, cur_val)
            cur_val += lazutkin_mesh
        return


    def import_fourier(self, path):
        """ Parses a text file located at path containing fourier coefficients 
            for the domain. Coefficients are expected one per line. 
        """

        self._clean_domain()
        file = open(path, "r")

        for line in file:
            self.fourier.append(mpf(line))
        arc_length = arc_length_coords(self, 1)
        self.fourier = [fdiv(coeff, arc_length) for coeff in self.fourier]
        
        curvature_minimum = minimize_scalar(
            lambda t: fdiv(1, self.radius_of_curv(t)), bracket=(0, 1 / 2)).fun
        curvature_maximum = fdiv(1, 
            minimize_scalar(lambda t: self.radius_of_curv(t), bracket=(0, 1 / 2)).fun)

        if (curvature_minimum < curvature_lower_bound):
            self.fourier = []
            raise Exception("The minimum curvature value should not be "
                             + "less than {}. The minimum value was {}."
                             .format(curvature_lower_bound, curvature_minimum))

        elif (curvature_maximum > curvature_upper_bound):
            self.fourier = []
            raise Exception("The maximum curvature value should not be "
                             + "larger than {}. The maximum value was {}."
                             .format(curvature_upper_bound, curvature_maximum))
        self._cache_lazutkin()
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

    def radius_second_derivative(self, theta):
        """ Computes the second derivative of the radius function at a given 
            angle. Assumed that import_fourier has been called. 
        """

        series_terms = [fprod([-1, power(fprod([2, pi, n]), 2), 
                               self.fourier[n]]) 
                        for n in range(len(self.fourier))]
        return fourierval((series_terms, [0]), [0, 1], theta)

    def radius_third_derivative(self, theta):
        """ Computes the second derivative of the radius function at a given
            angle. Assumed that import_fourier has been called
        """

        series_terms = [fmul(power(fprod([2, pi, n]), 3), 
                             self.fourier[n]) 
                        for n in range(len(self.fourier))]
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
        grad_x = fsub(fmul(r_prime, cos(theta_rescaled)), 
                      fprod([r, 2, pi, sin(theta_rescaled)]))
        grad_y = fadd(fmul(r_prime, sin(theta_rescaled)), 
                      fprod([r, 2, pi, cos(theta_rescaled)]))
        return (grad_x, grad_y)

    def polar_gradient_norm_deriv(self, theta):
        """ Returns the derivative of the norm of polar_gradient evaluated at 
            the point theta.
        """

        gradient = self.polar_gradient(theta)
        gradient_norm = norm(gradient)
        radius = self.radius(theta)
        radius_deriv = self.radius_derivative(theta)
        radius_second_deriv = self.radius_second_derivative(theta)
        theta_rescaled = fprod([theta, 2, pi])
        two_pi = fmul(2, pi)
        
        coeff_1 = fsub(radius_second_deriv, fmul(power(two_pi, 2), radius))
        coeff_2 = fprod([2, two_pi, radius_deriv])

        grad_x_prime = fsub(fmul(coeff_1, cos(theta_rescaled)), 
                            fmul(coeff_2, sin(theta_rescaled)))
        grad_y_prime = fadd(fmul(coeff_2, cos(theta_rescaled)), 
                            fmul(coeff_1, sin(theta_rescaled)))
        norm_deriv = fdiv(fadd(fmul(gradient[0], grad_x_prime), 
                               fmul(gradient[1], grad_y_prime)), 
                          gradient_norm)
        return norm_deriv

    def radius_of_curv(self, theta):
        """ Returns the radius of curvature at theta.
        """

        gradient = self.polar_gradient(theta)
        radius = self.radius(theta)
        radius_deriv = self.radius_derivative(theta)
        radius_second_deriv = self.radius_second_derivative(theta)
        theta_rescaled = fprod([theta, 2, pi])

        coeff_1 = fsub(radius_second_deriv, fprod([4, power(pi, 2), radius]))
        coeff_2 = fprod([4, pi, radius_deriv])

        second_grad_x = fsub(fmul(coeff_1, cos(theta_rescaled)), 
                             fmul(coeff_2, sin(theta_rescaled)))
        second_grad_y = fadd(fmul(coeff_2, cos(theta_rescaled)), 
                             fmul(coeff_1, sin(theta_rescaled)))
        denom = fsub(fmul(gradient[0], second_grad_y), 
                     fmul(gradient[1], second_grad_x))
        return abs(fdiv(power(norm(gradient), 3), denom))

    def radius_of_curv_deriv(self, theta):
        """ Returns the derivative of the radius of curvature at theta.
        """

        gradient = self.polar_gradient(theta)
        gradient_norm = norm(gradient)
        radius = self.radius(theta)
        radius_deriv = self.radius_derivative(theta)
        radius_second_deriv = self.radius_second_derivative(theta)
        radius_third_deriv = self.radius_third_derivative(theta)
        theta_rescaled = fprod([theta, 2, pi])

        coeff_1 = fsub(radius_second_deriv, fprod([4, power(pi, 2), radius]))
        coeff_2 = fprod([4, pi, radius_deriv])

        second_grad_x = fsub(fmul(coeff_1, cos(theta_rescaled)), 
                             fmul(coeff_2, sin(theta_rescaled)))
        second_grad_y = fadd(fmul(coeff_2, cos(theta_rescaled)), 
                             fmul(coeff_1, sin(theta_rescaled)))
        cross_product = fsub(fmul(gradient[0], second_grad_y), 
                             fmul(gradient[1], second_grad_x))

        coeff_1 = fsub(radius_third_deriv, fprod([12, power(pi, 2), 
                                                  radius_deriv]))
        coeff_2 = fmul(fmul(2, pi), fsub(fprod([4, power(pi, 2), radius]), 
                                         fmul(3, radius_second_deriv)))

        third_grad_x = fadd(fmul(coeff_1, cos(theta_rescaled)), 
                            fmul(coeff_2, sin(theta_rescaled)))
        third_grad_y = fsub(fmul(coeff_1, sin(theta_rescaled)),
                            fmul(coeff_2, cos(theta_rescaled)))

        numerator_term_1 = fprod([3, gradient_norm, 
                                  fadd(fmul(gradient[0], second_grad_x), 
                                       fmul(gradient[1], second_grad_y))])
        numerator_term_1 = fdiv(numerator_term_1, cross_product)

        numerator_term_2 = fprod([power(gradient_norm, 3), 
                                  fsub(fmul(gradient[0], third_grad_y), 
                                       fmul(gradient[1], third_grad_x))])
        numerator_term_2 = fdiv(numerator_term_2, power(cross_product, 2))

        numerator = fsub(numerator_term_1, numerator_term_2)
        denom = fdiv(cross_product, abs(cross_product))
        return fdiv(numerator, denom)

    def radius_of_curv_second_deriv(self, theta):
        """ Returns the second derivative of the radius of curvature at theta.
        """

        return diff(lambda t: self.radius_of_curv_deriv(t), theta, 
                    tol=power(10, -2 * mp.dps))