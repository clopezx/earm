"""Model from Chen 2007, FEBS Letters."""

from pysb import *
from earm2 import macros
from earm2 import earm2_modules
from earm2 import albeck_modules
from earm2 import shen_modules
from pysb.integrate import odesolve
from pysb.bng import generate_equations
from pylab import linspace, plot, figure, ion, legend
import re

Model()

shen_modules.momp_monomers()
#shen_modules.momp_initial_conditions(bid_state='T')

# The specific MOMP model to use
shen_modules.chen2007FEBS_indirect(pore_assembly=True)

Observable('Bax4_', MatchOnce(Bax(bf=None, s1=1, s2=4) %
                              Bax(bf=None, s1=2, s2=1) %
                              Bax(bf=None, s1=3, s2=2) %
                              Bax(bf=None, s1=4, s2=3)))
Observable('Bid_', Bid(bf=None))
Observable('Bcl2_', Bcl2(bf=None))
Observable('Bcl2_Bid_', Bcl2(bf=1) % Bid(bf=1))
Observable('Bcl2_Bax_', Bcl2(bf=1) % Bax(bf=1))

def figure_2a():
    f_range = linspace(0, 20, 50)
    t = linspace(0, 50000, 1000)
    ion()
    figure()
    ss_Bax4_vals = []

    for f in f_range:
        model.parameters['Bid_0'].value = 3 + (f * 3)
        x = odesolve(model, t)
        Bax_frac = (4*x['Bax4_'])/model.parameters['Bax_0'].value
        plot(t, Bax_frac, label='Bax4')
        ss_Bax4_val = Bax_frac[-1]
        ss_Bax4_vals.append(ss_Bax4_val)

    figure()
    plot(f_range, ss_Bax4_vals)


def print_odes():
    """Examine the structure of the equations after converting the nomenclature
       to the one used in the paper"""

    p_name_map = {
        'bind_BidT_Bcl2_kf': 'k_BH3_Bcl2',
        'bind_BidT_Bcl2_kr': 'kr_BH3Bcl2',
        'bind_BaxC_Bcl2_kf': 'k_Bax_Bcl2',
        'bind_BaxC_Bcl2_kr': 'kr_BaxBcl2',
        'spontaneous_pore_BaxC_to_Bax4_kf': 'k_o',
        'spontaneous_pore_BaxC_to_Bax4_kr': 'kr_o'
    }

    s_name_map = {'s0': 'BH3',
                  's1': 'Bax',
                  's2': 'Bcl2',
                  's3': 'BH3Bcl2',
                  's4': 'BaxBcl2',
                  's5': 'MAC'}

    generate_equations(model)

    for i, ode in enumerate(model.odes):
        new_ode = 'd[s%d]/dt = %s' % (i, str(ode))

        for old_s in s_name_map:
            new_ode = re.sub(old_s, s_name_map[old_s], new_ode)
        for old_p in p_name_map:
            new_ode = re.sub(old_p, p_name_map[old_p], new_ode)

        print new_ode




"""def simulate(f=0.5, tmax=3000):
    t = linspace(0, tmax, 100)  
    ion()
    figure()
    model.parameters['Bid_0'].value = f * 3
    x = odesolve(model, t)
    plot(t, (4*x['Bax4_'])/model.parameters['Bax_0'].value, label='Bax4')
    plot(t, x['Bid_']/model.parameters['Bid_0'].value, label='Bid free')
    plot(t, x['Bcl2_']/model.parameters['Bcl2_0'].value, label='Bcl2 free')
    plot(t, x['Bcl2_Bid_']/model.parameters['Bid_0'].value, label='Bcl2-Bid')
    plot(t, x['Bcl2_Bax_']/model.parameters['Bax_0'].value, label='Bcl2-Bax')
    legend()
    return x
"""



