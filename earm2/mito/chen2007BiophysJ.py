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
shen_modules.chen2007BiophysJ(pore_assembly=True)

Observable('aBax_', Bax(state='A', bf=None))
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
    f2_range = linspace(0, 1, 20)
    t = linspace(0, 50000, 1000)
    ion()
    figure()
    ss_aBax_vals = []

    for f2 in f2_range:
        model.parameters['Bid_0'].value = f2 * 1e-1
        x = odesolve(model, t)
        plot(t, x['aBax_']/model.parameters['Bax_0'].value)
        ss_aBax_val = x['aBax_'][-1]/model.parameters['Bax_0'].value
        ss_aBax_vals.append(ss_aBax_val)

    figure()
    plot(f2_range, ss_aBax_vals)

def print_original_odes():
    """Show the ODEs using the nomenclature from the paper.

    Example
    =======
    >>> from earm2.mito.chen2007BiophysJ import print_original_odes
    >>> print_original_odes()
    d[Act]/dt = -k5*Act*Bcl2 + k6*ActBcl2 + k7*AcBax*ActBcl2 - k8*Act*AcBaxBcl2
    d[InBax]/dt = -k1*Act*InBax + k2*AcBax
    d[Bcl2]/dt = -k3*Bcl2*AcBax + k4*AcBaxBcl2 - k5*Act*Bcl2 + k6*ActBcl2
    d[AcBax]/dt = -k3*Bcl2*AcBax + k4*AcBaxBcl2 - k7*AcBax*ActBcl2 + k8*Act*AcBaxBcl2 + k1*Act*InBax - k2*AcBax - 1.0*AcBax**4*k9 + 4*Bax4*k10
    d[ActBcl2]/dt = k5*Act*Bcl2 - k6*ActBcl2 - k7*AcBax*ActBcl2 + k8*Act*AcBaxBcl2
    d[AcBaxBcl2]/dt = k3*Bcl2*AcBax - k4*AcBaxBcl2 + k7*AcBax*ActBcl2 - k8*Act*AcBaxBcl2
    d[Bax4]/dt = 0.25*AcBax**4*k9 - Bax4*k10
    """

    generate_equations(model)

    # Examine the structure of the equations after converting the nomenclature
    # to the one used in the paper
    p_name_map = {
        'one_step_BidT_BaxC_to_BidT_BaxA_kf': 'k1',
        'reverse_BaxA_to_BaxC_k': 'k2',
        'bind_BidT_Bcl2_kf': 'k5',
        'bind_BidT_Bcl2_kr': 'k6',
        'bind_BaxA_Bcl2_kf': 'k3',
        'bind_BaxA_Bcl2_kr': 'k4',
        'displace_BaxA_BidTBcl2_to_BaxABcl2_BidT_kf': 'k7',
        'displace_BaxA_BidTBcl2_to_BaxABcl2_BidT_kr': 'k8',
        'spontaneous_pore_BaxA_to_Bax4_kf': 'k9',
        'spontaneous_pore_BaxA_to_Bax4_kr': 'k10'
    }

    s_name_map = {'s0': 'Act',
                  's1': 'InBax',
                  's2': 'Bcl2',
                  's3': 'AcBax',
                  's4': 'ActBcl2',
                  's5': 'AcBaxBcl2',
                  's6': 'Bax4'}

    for i, ode in enumerate(model.odes):
        new_ode = 'd[s%d]/dt = %s' % (i, str(ode))

        for old_s in s_name_map:
            new_ode = re.sub(old_s, s_name_map[old_s], new_ode)
        for old_p in p_name_map:
            new_ode = re.sub(old_p, p_name_map[old_p], new_ode)

        print new_ode

if __name__ == "__main__":
    import doctest
    doctest.testmod()
