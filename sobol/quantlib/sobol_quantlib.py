# Sample generation of 2 dimensions of Sobol sequence numbers using QuantLib
# source at: http://quantlib.org/
#
# I downloaded the source and the SWIG python files to use it as a module in python
#
#

import QuantLib
import numpy
import pylab

sobinit = QuantLib.SobolRsg(2, 64) # first number is number of dimensions (d), second is seed. They recommend seed of 64 (read in a forum)
                                   # this code can handle up to d< 360 dimensions

# now generate 1024 points. They recommend using powers of two for samples? (read in a forum)...
sobvals = [] # Hold the sobol values
for i in range(1024):
    t = sobinit.nextSequence() # returns a sequence object with a d-dimensional tuple (i.e. "d" values)
    sobvals.append(t.value())  # values from the sequence object

sobvals = numpy.asarray(sobvals)
sobvals = sobvals.T

pylab.plot(sobvals[0], sobvals[1], '.')
pylab.show()
