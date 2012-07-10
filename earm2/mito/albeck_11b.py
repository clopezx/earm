from pysb import *
from earm2 import albeck_modules
from sympy.parsing.sympy_parser import parse_expr

Model()

albeck_modules.momp_monomers()

# The specific MOMP model to use
albeck_modules.albeck_11b(do_pore_transport=True)

# Observables
Observable('AcBax_', Bax(bf=None, state='A'))
Observable('Bid_', Bid(bf=None))
#Observable('Bcl2_', Bcl2(bf=None))
#Observable('Bcl2_Bid_', Bcl2(bf=1) % Bid(bf=1))
Observable('Bcl2_Bax_', Bcl2(bf=1) % Bax(bf=1))
Observable('cSmac', Smac(state='C'))

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


