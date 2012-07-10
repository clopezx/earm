"""Model from Albeck et al., 2008, Figure 11b.
Minimal Model.
"""

from pysb import *
from earm2 import albeck_modules
from pysb.bng import generate_equations
from earm2.macros import catalyze
from pysb.integrate import odesolve
from matplotlib.pyplot import figure, ion, plot, legend
import numpy as np
from sympy.parsing.sympy_parser import parse_expr

rtol = 1e-6

#import albeck_11b
#import albeck_11c
import albeck_11e

m = albeck_11e.model
p_name_map = albeck_11e.p_name_map
s_index_map = albeck_11e.s_index_map
m_ode_list = albeck_11e.m_ode_list


def add_caspase8(model):
    # Add upstream caspase reaction to model
    if model.monomers.get('C8') is None:
        Bid = model.monomers.get('Bid')
        C8 = Monomer('C8', ['state', 'bf'], {'state': ['pro', 'A']},
                     _export=False)
        model.add_component(C8)
        C8_0 = Parameter('C8_0', 1, _export=False)
        model.add_component(C8_0)
        model.initial(C8(state='A', bf=None), C8_0)
        for component in catalyze(C8(state='A'), Bid(state='U'),
                                Bid(state='T'), [1e-7, 1e-3, 1]):
            model.add_component(component)

    # Set CytoC to 0 so transport is only of Smac
    #model.parameters['CytoC_0'].value = 0

def run_figure_sim(model):
    add_caspase8(model)

    # Set timepoints and c8 doses
    tf = 15 * 3600 # 15 hours
    t = np.linspace(0, tf, (tf / 60) + 1)
    c8_doses = [0.01e2, 0.05e2, 0.1e2, 0.5e2, 1e2, 5e2, 10e2]

    outputs = np.empty((len(t), len(c8_doses)))
    for i, c8_dose in enumerate(c8_doses):
        model.parameters['C8_0'].value = c8_dose
        x = odesolve(model, t, rtol=rtol)
        frac_Smac_release = x['cSmac'] / model.parameters['Smac_0'].value
        outputs[:,i] = frac_Smac_release

    return [t, outputs]

def plot_figure(model, data_file):
    [t, pysb_data] = run_figure_sim(model)
    figure()
    ion()
    pysb_num_doses = pysb_data.shape[1]

    mat_data = np.loadtxt(data_file)
    mat_num_doses = mat_data.shape[1]
    if pysb_num_doses != mat_num_doses:
        raise Exception('Different number of doses between PySB and MATLAB ' +
                        'files.')

    for i in range(pysb_num_doses):
        plot(t, pysb_data[:,i], 'r', label='PySB')
        plot(t, mat_data[:,i], 'b', label='Matlab')

    #legend(loc='lower right')


def matches_figure(model, data_file):
    [t, pysb_data] = run_figure_sim(model)
    mat_data = np.loadtxt(data_file)

    if pysb_data.shape[1] != mat_data.shape[1]:
        raise Exception('Different number of doses between PySB and MATLAB ' +
                        'files.')

    mat_data - pysb_data
    #import pdb; pdb.set_trace()
    #return bool((abs(diffs) < rtol * 10).all())
    return bool(np.allclose(pysb_data, mat_data, atol=rtol*10))


def compare_odes(model, p_name_map, s_index_map, m_ode_list):
    add_caspase8(model)
    s_name_map = [('s%d' % i, 'x(%d)' % j) for i, j in enumerate(s_index_map)]
    name_map = {}
    name_map.update(p_name_map)
    name_map.update(s_name_map)
    generate_equations(model)
    ode_list = []
    result_list = []
    for i, ode in enumerate(model.odes):
        new_ode = ode.subs(name_map)
        #new_ode = 'd[%s]/dt = %s' % (s_name_map['s%d' % i], str(new_ode))
        #ode_list.append(new_ode)
        old_ode = m_ode_list[s_index_map[i] - 1]
        result = new_ode == old_ode
        result_list.append(result)
        if not result:
            print "Mismatch for species " + str(i) + ": "
            print "Pysb ODE   : %s" % str(new_ode)
            print "Matlab ODE : %s" % str(old_ode)

    return result_list

compare_odes(m, p_name_map, s_index_map, m_ode_list)
plot_figure(m, 'albeck_11e.tsv')
print matches_figure(m, 'albeck_11e.tsv')
