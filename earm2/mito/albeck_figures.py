"""Model from Albeck et al., 2008, Figure 11b.
Minimal Model.
"""

from pysb import *
from earm2 import albeck_modules
from pysb.bng import generate_equations
from earm2.macros import catalyze
from pysb.integrate import odesolve
from matplotlib.pyplot import figure, ion, plot
import numpy as np
from sympy.parsing.sympy_parser import parse_expr

import albeck_11b

rtol = 1e-6
m = albeck_11b.model

def run_figure_sim(model):
    # Add upstream caspase reaction to model
    if model.monomers.get('C8') is None:
        Bid = model.monomers.get('Bid')
        C8 = Monomer('C8', ['state', 'bf'], {'state': ['pro', 'A']})
        model.add_component(C8)
        C8_0 = Parameter('C8_0', 0)
        model.add_component(C8_0)
        model.initial(C8(state='A', bf=None), C8_0)
        for component in catalyze(C8(state='A'), Bid(state='U'),
                                Bid(state='T'), [1e-7, 1e-3, 1]):
            model.add_component(component)

    # Set CytoC to 0 so transport is only of Smac
    model.parameters['CytoC_0'].value = 0

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
        plot(t, pysb_data[:,i], 'r')
        plot(t, mat_data[:,i], 'b')

def matches_figure(model, data_file):
    [t, pysb_data] = run_figure_sim(model)
    mat_data = np.loadtxt(data_file)

    if pysb_data.shape[1] != mat_data.shape[1]:
        raise Exception('Different number of doses between PySB and MATLAB ' +
                        'files.')

    diffs = mat_data - pysb_data
    return (diffs < rtol * 10).all()

def convert_odes():
    p_name_map = {
        'bind_BidT_BaxC_to_BidTBaxC_kf': 'k(2)',
        'bind_BidT_BaxC_to_BidTBaxC_kr': 'k_(2)',
        'catalyze_BidTBaxC_to_BidT_BaxA_kc': 'kc(2)',
        'bind_BidT_Bcl2_kf': 'k(4)',
        'bind_BidT_Bcl2_kr': 'k_(4)',
        'bind_BaxA_Bcl2_kf': 'k(5)',
        'bind_BaxA_Bcl2_kr': 'k_(5)',
        'bind_BaxA_SmacM_to_BaxASmacM_kf': 'k(3)',
        'bind_BaxA_SmacM_to_BaxASmacM_kr': 'k_(3)',
        'catalyze_BaxASmacM_to_BaxA_SmacC_kc': 'kc(3)',
        'bind_C8A_BidU_to_C8ABidU_kf': 'k(1)',
        'bind_C8A_BidU_to_C8ABidU_kr': 'k_(1)',
        'catalyze_C8ABidU_to_C8A_BidT_kc': 'kc(1)'}
    s_index_map = [2, 4, 8, 6, 1, 9, 3, 10, 12, 5, 13, 11, 7]
    s_name_map = [('s%d' % i, 'x(%d)' % j) for i, j in enumerate(s_index_map)]

    # The list of ODEs from John Albeck's Matlab file
    m_ode_list = [
        '-k(1)*x(1)*x(2) +k_(1)*x(9)+kc(1)*x(9)',
        '-k(1)*x(1)*x(2) +k_(1)*x(9)',
        'kc(1)*x(9) -k(2)*x(3)*x(4) +k_(2)*x(10)+kc(2)*x(10) - k(4)*x(3)*x(8)+ k_(5)*x(12)',
        '-k(2)*x(3)*x(4)+k_(2)*x(10)',
        'kc(2)*x(10) -  k(3)*x(5)*x(6)  + k_(3)*x(11) +kc(3)*x(11) - k(5)*x(5)*x(8)+ k_(5)*x(13)',
        '-k(3)*x(5)*x(6) +k_(3)*x(11)',
        'kc(3)*x(11)',
        '-k(4)*x(3)*x(8)+k_(4)*x(12) - k(5)*x(5)*x(8)+k_(5)*x(13)',
        'k(1)*x(1)*x(2)-k_(1)*x(9)-kc(1)*x(9)',
        'k(2)*x(3)*x(4)-k_(2)*x(10)-kc(2)*x(10)',
        'k(3)*x(5)*x(6)-k_(3)*x(11)-kc(3)*x(11)',
        'k(4)*x(3)*x(8)-k_(4)*x(12)',
        'k(5)*x(5)*x(8)-k_(5)*x(13)']
    m_ode_list = [parse_expr(ode) for ode in m_ode_list]

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
        if not result:
            print new_ode
            print old_ode

    return result_list

