"""Model from Cui 2008, PLoS One.
TODO: Docstring
"""

from pysb import *
from earm2 import shen_modules
from pysb.bng import generate_equations
import re

Model()

shen_modules.momp_monomers()

# The specific MOMP model to use
shen_modules.cui2008_direct()

# Observables
Observable('AcBax_', Bax(bf=None, state='A'))
#Observable('Bax4_', MatchOnce(Bax(bf=None, s1=1, s2=4) %
#                              Bax(bf=None, s1=2, s2=1) %
#                              Bax(bf=None, s1=3, s2=2) %
#                              Bax(bf=None, s1=4, s2=3)))
Observable('Bid_', Bid(bf=None))
Observable('Bad_', Bad(bf=None))
Observable('Bcl2_', Bcl2(bf=None))
Observable('Bcl2_Bad_', Bcl2(bf=1) % Bad(bf=1))
Observable('Bcl2_Bid_', Bcl2(bf=1) % Bid(bf=1))
Observable('Bcl2_Bax_', Bcl2(bf=1) % Bax(bf=1))

def print_original_odes():
    """Show the ODEs using the nomenclature from the paper.

    Example
    -------
    Note that the coefficient of 0.5 in the final equation
    (for the pore, a.k.a. MAC) that is introduced by BNG is accounted
    for by multiplying the forward rate of pore formation, k16,
    by 2 (see cui2008_direct() in shen_modules.py):

    >>> import earm2.mito.cui2008_direct as c
    >>> c.print_original_odes()
    d[Act]/dt = -k4*Act*Bcl2 + k5*ActBcl2 - u3*Act + k11*Ena*ActBcl2 - k12*Act*EnaBcl2 + k6*AcBax*ActBcl2 - k7*Act*AcBaxBcl2 + __source*p2
    d[Ena]/dt = -k9*Ena*Bcl2 + k10*EnaBcl2 - u7*Ena - k13*Ena*AcBaxBcl2 + k14*AcBax*EnaBcl2 - k11*Ena*ActBcl2 + k12*Act*EnaBcl2 + __source*p4
    d[InBax]/dt = -u1*InBax - k1*Act*InBax + k8*AcBax + __source*p1
    d[Bcl2]/dt = -k9*Ena*Bcl2 + k10*EnaBcl2 - k4*Act*Bcl2 + k5*ActBcl2 - u4*Bcl2 + __source*p3
    d[__source]/dt = 0
    d[AcBax]/dt = -u2*AcBax - 1.0*k16*AcBax**2 + 2*k17*MAC + k13*Ena*AcBaxBcl2 - k14*AcBax*EnaBcl2 - k6*AcBax*ActBcl2 + k7*Act*AcBaxBcl2 + k1*Act*InBax - k8*AcBax
    d[ActBcl2]/dt = k4*Act*Bcl2 - k5*ActBcl2 - u5*ActBcl2 - k11*Ena*ActBcl2 + k12*Act*EnaBcl2 - k6*AcBax*ActBcl2 + k7*Act*AcBaxBcl2
    d[EnaBcl2]/dt = k9*Ena*Bcl2 - k10*EnaBcl2 - u8*EnaBcl2 + k13*Ena*AcBaxBcl2 - k14*AcBax*EnaBcl2 + k11*Ena*ActBcl2 - k12*Act*EnaBcl2
    d[__sink]/dt = u8*EnaBcl2 + u7*Ena + u2*AcBax + u9*MAC + u6*AcBaxBcl2 + u1*InBax + u4*Bcl2 + u5*ActBcl2 + u3*Act
    d[AcBaxBcl2]/dt = -u6*AcBaxBcl2 - k13*Ena*AcBaxBcl2 + k14*AcBax*EnaBcl2 + k6*AcBax*ActBcl2 - k7*Act*AcBaxBcl2
    d[MAC]/dt = -u9*MAC + 0.5*k16*AcBax**2 - k17*MAC
    """

    p_name_map = {
        'one_step_BidT_BaxC_to_BidT_BaxA_kf': 'k1',
        'reverse_BaxA_to_BaxC_k': 'k8',
        'bind_BidT_Bcl2_kf': 'k4',
        'bind_BidT_Bcl2_kr': 'k5',
        'bind_BadM_Bcl2_kf': 'k9',
        'bind_BadM_Bcl2_kr': 'k10',
        'displace_BaxA_BidTBcl2_to_BaxABcl2_BidT_kf': 'k6',
        'displace_BaxA_BidTBcl2_to_BaxABcl2_BidT_kr': 'k7',
        'displace_BadM_BidTBcl2_to_BadMBcl2_BidT_kf': 'k11',
        'displace_BadM_BidTBcl2_to_BadMBcl2_BidT_kr': 'k12',
        'displace_BadM_BaxABcl2_to_BadMBcl2_BaxA_kf': 'k13',
        'displace_BadM_BaxABcl2_to_BadMBcl2_BaxA_kr': 'k14',
        'dimerize_Bax_kf': 'k16',
        'dimerize_Bax_kr': 'k17',
        'synthesize_BaxC_k': 'p1',
        'degrade_BaxC_k': 'u1',
        'degrade_BaxA_k': 'u2',
        'synthesize_BidT_k': 'p2',
        'degrade_BidT_k': 'u3',
        'synthesize_Bcl2_k': 'p3',
        'degrade_Bcl2_k': 'u4',
        'degrade_BidTBcl2_k': 'u5',
        'degrade_BaxBcl2_k': 'u6',
        'synthesize_BadM_k': 'p4',
        'degrade_BadM_k': 'u7',
        'degrade_BadBcl2_k': 'u8',
        'degrade_BaxBax_k': 'u9'
    }

    s_name_map = {'s0': 'Act',
                  's1': 'Ena',
                  's2': 'InBax',
                  's3': 'Bcl2',
                  's4': '__source',
                  's5': 'AcBax',
                  's6': 'ActBcl2',
                  's7': 'EnaBcl2',
                  's8': '__sink',
                  's9': 'AcBaxBcl2',
                  's10': 'MAC'}

    generate_equations(model)

    for i, ode in enumerate(model.odes):
        new_ode = 'd[s%d]/dt = %s' % (i, str(ode))

        for old_s in s_name_map:
            new_ode = re.sub(old_s, s_name_map[old_s], new_ode)
        for old_p in p_name_map:
            new_ode = re.sub(old_p, p_name_map[old_p], new_ode)

        #if not (i == 4 or i == 8): # ignore __sink ODE
        print new_ode


if __name__ == "__main__":
    import doctest
    doctest.testmod()
