from pysb import *
from earm2 import albeck_modules
from sympy.parsing.sympy_parser import parse_expr

Model()

albeck_modules.momp_monomers()

# The specific MOMP model to use
albeck_modules.albeck_11d(do_pore_transport=True)

# Observables
Observable('AcBax_', Bax(bf=None, state='A'))
Observable('Bid_', Bid(bf=None))
#Observable('Bcl2_', Bcl2(bf=None))
#Observable('Bcl2_Bid_', Bcl2(bf=1) % Bid(bf=1))
Observable('Bcl2_Bax_', Bcl2(bf=1) % Bax(bf=1))
Observable('cSmac', Smac(state='C'))

p_name_map = {
    'bind_C8A_BidU_to_C8ABidU_kf': 'k(1)',
    'bind_C8A_BidU_to_C8ABidU_kr': 'k_(1)',
    'catalyze_C8ABidU_to_C8A_BidT_kc': 'kc(1)',
    'bind_BidT_BaxC_to_BidTBaxC_kf': 'k(2)',
    'bind_BidT_BaxC_to_BidTBaxC_kr': 'k_(2)',
    'catalyze_BidTBaxC_to_BidT_BaxA_kc': 'kc(2)',
    'equilibrate_BaxA_to_BaxM_kf': 'k(3)',
    'equilibrate_BaxA_to_BaxM_kr': 'k_(3)',
    'Bax_dimerization_kf': 'k(4)',
    'Bax_dimerization_kr': 'k_(4)',
    'Bax_tetramerization_kf': 'k(5)',
    'Bax_tetramerization_kr': 'k_(5)',
    'pore_bind_Bax_2_Bcl2_kf': 'k(8)',
    'pore_bind_Bax_2_Bcl2_kr': 'k_(8)',
    'pore_bind_Bax_4_Bcl2_kf': 'k(9)',
    'pore_bind_Bax_4_Bcl2_kr': 'k_(9)',
    'bind_BaxM_Bcl2_kf': 'k(7)',
    'bind_BaxM_Bcl2_kr': 'k_(7)',
    'pore_transport_complex_Bax_4_SmacM_kf': 'k(6)',
    'pore_transport_complex_Bax_4_SmacM_kr': 'k_(6)',
    'pore_transport_dissociate_Bax_4_SmacC_kc': 'kc(6)'}

s_index_map = [2, 4, 11, 9, 1, 12, 3, 13, 5, 6, 7, 15, 8, 16, 17, 14, 10]


# The list of ODEs from John Albeck's Matlab file
m_ode_list = [
    '-k(1)*x(1)*x(2) +k_(1)*x(12)+kc(1)*x(12)', # cCasp8
    '-k(1)*x(1)*x(2) +k_(1)*x(12)', # Bid
    'kc(1)*x(12) -k(2)*x(3)*x(4) +k_(2)*x(13)+kc(2)*x(13)' # tBid
    '-k(2)*x(3)*x(4)+k_(2)*x(13)', # Bax
    'kc(2)*x(13)-k(3)*x(5)+k_(3)*x(6)', # aBax
    'k(3)*x(5)-k_(3)*x(6) -1/v*2*k(4)*x(6)**2+2*k_(4)*x(7) -1/v*k(7)*x(6)*x(11)+k_(7)*x(15)', # mBax
    '1/v*k(4)*x(6)**2-k_(4)*x(7) -1/v*2*k(5)*x(7)**2+2*k_(5)*x(8) -1/v*k(8)*x(7)*x(11)+k_(8)*x(16)', # mBax2
    '1/v*k(5)*x(7)**2 - k_(5)*x(8) -1/v*k(6)*x(8)*x(9) +k_(6)*x(14)+kc(6)*x(14) -1/v*k(9)*x(8)*x(11)+k_(9)*x(17)', # mBax4
    '- 1/v*k(6)*x(8)*x(9) +k_(6)*x(14)', # Smac
    'kc(6)*x(14)', # aSmac
    '-1/v*k(7)*x(6)*x(11)+k_(7)*x(15) -1/v*k(8)*x(7)*x(11)+k_(8)*x(16) -1/v*k(9)*x(8)*x(11)+k_(9)*x(17)', # Bcl2
    'k(1)*x(1)*x(2)-k_(1)*x(12)-kc(1)*x(12)', # cCasp8:Bid
    'k(2)*x(3)*x(4)-k_(2)*x(13)-kc(2)*x(13)', # tBid:Bax
    '1/v*k(6)*x(8)*x(9)-k_(6)*x(14)-kc(6)*x(14)', # mBax4:Smac
    '1/v*k(7)*x(6)*x(11)-k_(7)*x(15)', # aBax:Bcl2
    '1/v*k(8)*x(7)*x(11)-k_(8)*x(16)', # aBax:Bcl2
    '1/v*k(9)*x(8)*x(11)-k_(9)*x(17)', # aBax:Bcl2
]
m_ode_list = [parse_expr(ode) for ode in m_ode_list]
