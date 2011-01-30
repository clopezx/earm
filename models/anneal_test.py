import scipy.optimize
import numpy

# From the Landau spinodal theory:
# Define
# V(x) = (1/2)rx^2 + (1/4)ux^4
# Then for r>0 you have basically a parabola with minimum at x=0
# and for r<0 you have double well potential, where x=0 is the local maximum (the
# barrier between the wells) and the well locations are at + or - sqrt(|r|/u)
def double_well(x, r = -7., u = 1.):
    return 0.5*(r-0.2*x)*(x**2) + 0.25*u*(x**4)
    

xdata = numpy.arange(-4, 5, 0.1)

for i in xdata:
    ydata = double_well(xdata)

