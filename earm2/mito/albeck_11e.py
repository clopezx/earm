"""TODO: Docstring"""
from pysb import *
from earm2 import albeck_modules
from sympy.parsing.sympy_parser import parse_expr

Model()

albeck_modules.momp_monomers()

# The specific MOMP model to use
albeck_modules.albeck_11e(do_pore_transport=True)

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
    'equilibrate_BaxA_to_BaxM_kf': 'k(3)',
    'equilibrate_BaxA_to_BaxM_kr': 'k_(3)',
    'Bax_dimerization_kf': 'k(4)',
    'Bax_dimerization_kr': 'k_(4)',
    'Bax_tetramerization_kf': 'k(5)',
    'Bax_tetramerization_kr': 'k_(5)',
    'pore_bind_Bax_2_Bcl2_kf': 'k(9)',
    'pore_bind_Bax_2_Bcl2_kr': 'k_(9)',
    'pore_bind_Bax_4_Bcl2_kf': 'k(10)',
    'pore_bind_Bax_4_Bcl2_kr': 'k_(10)',
    'bind_BaxM_Bcl2_kf': 'k(8)',
    'bind_BaxM_Bcl2_kr': 'k_(8)',
    'bind_MitoA_SmacM_to_MitoASmacM_kf': 'k(7)',
    'bind_MitoA_SmacM_to_MitoASmacM_kr': 'k_(7)',
    'catalyze_MitoASmacM_to_MitoA_SmacC_kc': 'kc(7)',
    'pore_bind_Bax_4_MitoI_kf': 'k(6)',
    'pore_bind_Bax_4_MitoI_kr': 'k_(6)',
    'Mito_activation_kc': 'kc(6)'}
s_index_map = [2, 4, 13, 9, 11, 1, 14, 3, 15, 5, 6, 7, 18, 8, 19, 20, 16, 10,
               17, 12]

# The list of ODEs from John Albeck's Matlab file
m_ode_list = [
    '-k(1)*x(1)*x(2) +k_(1)*x(14)+kc(1)*x(14)', # cCasp8
    '-k(1)*x(1)*x(2) +k_(1)*x(14)', # Bid
    'kc(1)*x(14) -k(2)*x(3)*x(4) +k_(2)*x(15)+kc(2)*x(15)', # tBid
    '-k(2)*x(3)*x(4)+k_(2)*x(15)', # Bax
    'kc(2)*x(15)-k(3)*x(5)+k_(3)*x(6) ', # aBax
    'k(3)*x(5)-k_(3)*x(6) -1/v*2*k(4)*x(6)**2+2*k_(4)*x(7) -1/v*k(8)*x(6)*x(13)+k_(8)*x(18)', # mBax
    '1/v*k(4)*x(6)**2-k_(4)*x(7) -1/v*2*k(5)*x(7)**2+2*k_(5)*x(8) -1/v*k(9)*x(7)*x(13)+k_(9)*x(19)', # mBax2
    '1/v*k(5)*x(7)**2 - k_(5)*x(8) - 1/v*k(6)*x(8)*x(9) +k_(6)*x(16) -1/v*k(10)*x(8)*x(13)+k_(10)*x(20)', # mBax4
    '- 1/v*k(6)*x(8)*x(9) +k_(6)*x(16)', # M
    'kc(6)*x(16) - 1/v*k(7)*x(10)*x(11)+k_(7)*x(17)+kc(7)*x(17)', # M*
    '- 1/v*k(7)*x(10)*x(11) +k_(7)*x(17)', # Smac
    'kc(7)*x(17)', # aSmac
    '-1/v*k(8)*x(6)*x(13) +k_(8)*x(18) -1/v*k(9)*x(7)*x(13) +k_(9)*x(19) -1/v*k(10)*x(8)*x(13)+k_(10)*x(20)', # Bcl2
    'k(1)*x(1)*x(2)-k_(1)*x(14)-kc(1)*x(14)', # cCasp8:Bid
    'k(2)*x(3)*x(4)-k_(2)*x(15)-kc(2)*x(15)', # tBid:Bax
    '1/v*k(6)*x(8)*x(9)-k_(6)*x(16)-kc(6)*x(16)', # mBax4:M
    '1/v*k(7)*x(10)*x(11) - k_(7)*x(17) - kc(7)*x(17)', # M*:Smac
    '1/v*k(8)*x(6)*x(13)-k_(8)*x(18)', # aBax:Bcl2
    '1/v*k(9)*x(7)*x(13)-k_(9)*x(19)', # aBax2:Bcl2
    '1/v*k(10)*x(8)*x(13)-k_(10)*x(20)', # aBax4:Bcl2
]
m_ode_list = [parse_expr(ode) for ode in m_ode_list]
