from pysb import *
from pysb.integrate import odesolve
from pylab import plot, linspace, xlabel, ylabel
import numpy.random

Model()

# Declare the monomers
Monomer('L', ['s'])
Monomer('R', ['s'])

# Declare the parameters
Parameter('kf', 1e-3)
Parameter('kr', 1e-3)

# Declare the initial conditions
Initial(L(s=None), Parameter('L_0', 100))
Initial(R(s=None), Parameter('R_0', 200))

# Declare the binding rule
Rule('L_binds_R',     
     L(s=None) + R(s=None) <>
     L(s=1) % R(s=1), kf, kr)

# Observe the complex
Observable('LR', L(s=1) % R(s=1))

# Simulate the model
time = linspace(0, 40, 100)
x = odesolve(model, time)
plot(time, x['LR'])

xlabel('Time (seconds)')
ylabel('Amount of LR')

def synthetic_data(tspan):
    seed = 1
    random = numpy.random.RandomState(seed)
    sigma = 0.1
    ntimes = 20
    #tspan = linspace(0, 40, ntimes)
    x = odesolve(model, tspan)
    lr_data = x['LR'] 
    lr_data = lr_data * (random.randn(len(lr_data)) * sigma + 1);
    #ysim_array = ysim.view().reshape(len(tspan), len(ysim.dtype))
    #yspecies = ysim_array[:, :len(model.species)]
    #ydata = yspecies * (random.randn(*yspecies.shape) * sigma + 1);
    #ysim_max = yspecies.max(0)
    #ydata_norm = ydata / ysim_max
    return lr_data
