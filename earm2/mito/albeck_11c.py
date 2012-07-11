"""TODO: Docstring"""
from pysb import *
from earm2 import albeck_modules
from sympy.parsing.sympy_parser import parse_expr

Model()

albeck_modules.momp_monomers()

# The specific MOMP model to use
albeck_modules.albeck_11c(do_pore_transport=True)

# Observables
Observable('AcBax_', Bax(bf=None, state='A'))
Observable('cSmac', Smac(state='C'))

p_name_map = {
    'bind_C8A_BidU_to_C8ABidU_kf': 'k(1)',
    'bind_C8A_BidU_to_C8ABidU_kr': 'k_(1)',
    'catalyze_C8ABidU_to_C8A_BidT_kc': 'kc(1)',
    'bind_BidT_BaxC_to_BidTBaxC_kf': 'k(2)',
    'bind_BidT_BaxC_to_BidTBaxC_kr': 'k_(2)',
    'catalyze_BidTBaxC_to_BidT_BaxA_kc': 'kc(2)',
    'Bax_dimerization_kf': 'k(3)',
    'Bax_dimerization_kr': 'k_(3)',
    'Bax_tetramerization_kf': 'k(4)',
    'Bax_tetramerization_kr': 'k_(4)',
    'pore_bind_Bax_2_Bcl2_kf': 'k(7)',
    'pore_bind_Bax_2_Bcl2_kr': 'k_(7)',
    'pore_bind_Bax_4_Bcl2_kf': 'k(8)',
    'pore_bind_Bax_4_Bcl2_kr': 'k_(8)',
    'bind_BaxA_Bcl2_kf': 'k(6)',
    'bind_BaxA_Bcl2_kr': 'k_(6)',
    'pore_transport_complex_Bax_4_SmacM_kf': 'k(5)',
    'pore_transport_complex_Bax_4_SmacM_kr': 'k_(5)',
    'pore_transport_dissociate_Bax_4_SmacC_kc': 'kc(5)'}

s_index_map = [2, 4, 10, 8, 1, 11, 3, 12, 5, 6, 14, 7, 15, 16, 13, 9]

# The list of ODEs from John Albeck's Matlab file
m_ode_list = [
    '-k(1)*x(1)*x(2) +k_(1)*x(11)+kc(1)*x(11)', # cCasp8
    '-k(1)*x(1)*x(2) +k_(1)*x(11)', # Bid
    'kc(1)*x(11) -k(2)*x(3)*x(4) +k_(2)*x(12)+kc(2)*x(12)', # tBid
    '-k(2)*x(3)*x(4)+k_(2)*x(12)', # Bax
    'kc(2)*x(12)- 2*k(3)*x(5)**2+2*k_(3)*x(6)- k(6)*x(5)*x(10)+ k_(6)*x(14)', # aBax
    'k(3)*x(5)**2-k_(3)*x(6)- 2*k(4)*x(6)**2+2*k_(4)*x(7)- k(7)*x(6)*x(10)+ k_(7)*x(15)', # aBax2
    'k(4)*x(6)**2 - k_(4)*x(7) -  k(5)*x(7)*x(8)  +k_(5)*x(13) +kc(5)*x(13) - k(8)*x(7)*x(10)+ k_(8)*x(16)', # aBax4
    '-k(5)*x(7)*x(8) +k_(5)*x(13)', # Smac
    'kc(5)*x(13)', # aSmac
    '-k(6)*x(5)*x(10)+k_(6)*x(14) - k(7)*x(6)*x(10)+k_(7)*x(15)- k(8)*x(7)*x(10)+k_(8)*x(16)', # Bcl2
    'k(1)*x(1)*x(2)-k_(1)*x(11)-kc(1)*x(11)', # cCasp8:Bid
    'k(2)*x(3)*x(4)-k_(2)*x(12)-kc(2)*x(12)',# tBid:Bax
    'k(5)*x(7)*x(8)-k_(5)*x(13)-kc(5)*x(13)', # aBax4:Smac
    'k(6)*x(5)*x(10)-k_(6)*x(14)', # aBax:Bcl2
    'k(7)*x(6)*x(10)-k_(7)*x(15)', # aBax2:Bcl2
    'k(8)*x(7)*x(10)-k_(8)*x(16)'] # aBax4:Bcl2
m_ode_list = [parse_expr(ode) for ode in m_ode_list]

