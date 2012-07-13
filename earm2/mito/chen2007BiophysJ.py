"""Model from Chen 2007, Biophys J."""

from pysb import *
from earm2 import macros
from earm2 import earm2_modules
from earm2 import albeck_modules
from earm2 import shen_modules
from pysb.integrate import odesolve
from pylab import linspace, plot, figure, ion, legend
from pysb.bng import generate_equations
import re

Model()

shen_modules.momp_monomers()
#shen_modules.momp_initial_conditions(bid_state='T')

# The specific MOMP model to use
shen_modules.chen2007BiophysJ(do_pore_assembly=True)

Observable('aBax_', Bax(state='A', bf=None))
#Observable('aBax_', Bax(state='A', bf=None, s1=None, s2=None))
#Observable('aBax_', Bax(state='A'))
Observable('Bid_', Bid(bf=None))
Observable('Bcl2_', Bcl2(bf=None))
Observable('Bcl2_Bid_', Bcl2(bf=1) % Bid(bf=1))
Observable('Bcl2_Bax_', Bcl2(bf=1) % Bax(bf=1))

def simulate(f2=0.5, tmax=3000):
    t = linspace(0, tmax, 100)  
    ion()
    figure()
    model.parameters['Bid_0'].value = f2 * 1e-1
    x = odesolve(model, t)
    plot(t, x['aBax_']/model.parameters['Bax_0'].value, label='aBax')
    plot(t, x['Bid_']/model.parameters['Bid_0'].value, label='Bid free')
    plot(t, x['Bcl2_']/model.parameters['Bcl2_0'].value, label='Bcl2 free')
    plot(t, x['Bcl2_Bid_']/model.parameters['Bid_0'].value, label='Bcl2-Bid')
    plot(t, x['Bcl2_Bax_']/model.parameters['Bax_0'].value, label='Bcl2-Bax')
    legend()
    return x

def bistability_analysis():
    f2_range = linspace(0, 0.4, 41)
    t = linspace(0, 50000, 1000)
    ion()
    ss_aBax_vals_up = []
    ss_aBax_vals_down = []

    for f2 in f2_range:
        model.parameters['Bid_0'].value = f2 * 1e-1
        bax_total = 2e-1

        # Do "up" portion of hysteresis plot
        model.parameters['aBax_0'].value = 0
        model.parameters['cBax_0'].value = bax_total
        x = odesolve(model, t)
        figure('up')
        plot(t, x['aBax_']/bax_total)
        ss_aBax_vals_up.append(x['aBax_'][-1]/bax_total)

        # Do "down" portion of hysteresis plot
        model.parameters['aBax_0'].value = bax_total
        model.parameters['cBax_0'].value = 0
        x = odesolve(model, t)
        figure('down')
        plot(t, x['aBax_']/bax_total)
        ss_aBax_vals_down.append(x['aBax_'][-1]/bax_total)

    figure()
    plot(f2_range, ss_aBax_vals_up, 'r')
    plot(f2_range, ss_aBax_vals_down, 'g')

