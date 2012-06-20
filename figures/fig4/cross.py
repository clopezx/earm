from pysb import *
from pysb import kappa
from pysb import bng
from pysb.integrate import odesolve
from pylab import figure, plot, show, linspace, xlabel, \
                  ylabel, title, xlim, ylim
from sympy import solve_poly_system, var, latex, S
from sympy.parsing.sympy_parser import parse_expr

from fig1.hello_pysb import model

# Instantiate the model
#Model()

#Monomer('tBid', ['b'])
#Monomer('Bax', ['b', 'state'], {'state': ['I','A']})

#two_step_mod(tBid(b=None), Bax(b=None, state='I'),
#             Bax(b=None, state='A'), site='b')

#Rule('Bax_inactivation',
#      Bax(b=None, state='A') >> Bax(b=None, state='I'),
#       Parameter('k_rev', 1e-1))

# Set the initial conditions
#Initial(tBid(b=None), Parameter('tBid_0', 1000))
#Initial(Bax(b=None, state='I'), Parameter('Bax_0', 1000))

# Model output: activated Bax
#Observe('aBax', Bax(state='A'))

# Run a deterministic ODE simulation using BNG and VODE
tmax = 40
numpoints = 100
det_time = linspace(0, tmax, numpoints)
det_data = odesolve(model, det_time)

#figure()
#plot(det_time, det_data['LR'])
#show()
"""
# Run a stochastic simulation using kappa
kappa_stoch_data = kappa.run_simulation(model, time=tmax, points=numpoints)

figure()
plot(kappa_stoch_data['time'], kappa_stoch_data['LR'])
xlim([0, 40])
ylim([0, 100])
title('Kappa Stochastic Simulation')
xlabel('Time')
ylabel('Amount of LR')
show()
                                        
# Run a stochastic simulation using BNG
bng_stoch_data = bng.run_ssa(model, t_end=tmax, n_steps=numpoints)

figure()
plot(bng_stoch_data['time'], bng_stoch_data['LR'])
xlim([0, 40])
ylim([0, 100])
title('BNG Stochastic Simulation')
xlabel('Time')
ylabel('Amount of LR')
show()

# Show contact and influence maps
#kappa.show_influence_map(model)
#kappa.show_contact_map(model)
"""

# Solve the model for steady state activated Bax using SymPy
#cons_eqns = parse_expr('s0 + s2 - tBid_0, 
bng.generate_equations(model)
var('s0, s1, s2, L_0, R_0')
conservation_eqns = [s0 + s2 - L_0, s1 + s2 - R_0]
solution = solve_poly_system(model.odes + conservation_eqns, s0, s1, s2)

s2_soln = solution[0][2].subs({S('kf'): S('k_f'), S('kr'): S('k_r')})
latex_output = latex(s2_soln)
print latex_output
