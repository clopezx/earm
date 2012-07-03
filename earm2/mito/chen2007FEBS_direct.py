"""Model from Chen 2007, FEBS Letters.
TODO: Docstring
"""

from pysb import *
from earm2 import shen_modules
from pysb.bng import generate_equations
import re

Model()

shen_modules.momp_monomers()

# The specific MOMP model to use
shen_modules.chen2007FEBS_direct(pore_assembly=True)

# Observables
Observable('AcBax_', Bax(bf=None, state='A'))
Observable('Bax4_', MatchOnce(Bax(bf=None, s1=1, s2=4) %
                              Bax(bf=None, s1=2, s2=1) %
                              Bax(bf=None, s1=3, s2=2) %
                              Bax(bf=None, s1=4, s2=3)))
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
    Note that the coefficient of 0.25 in the final equation
    (for the pore, a.k.a. MAC) that is introduced by BNG is accounted
    for by multiplying the forward rate of pore formation, k_o,
    by 4 (see chen2007FEBS_direct() in shen_modules.py):

    >>> import earm2.mito.chen2007FEBS_direct as d
    >>> d.print_original_odes()
    d[Act]/dt = -k_BH3_Bcl2*Act*Bcl2 + kr_BH3Bcl2*ActBcl2
    d[Ena]/dt = -k_BH3_Bcl2*Ena*Bcl2 + kr_BH3Bcl2*EnaBcl2
    d[InBax]/dt = -k_InBax*Act*InBax + k_Bax*Bax
    d[Bcl2]/dt = -k_BH3_Bcl2*Ena*Bcl2 + kr_BH3Bcl2*EnaBcl2 - k_BH3_Bcl2*Act*Bcl2 + kr_BH3Bcl2*ActBcl2
    d[Bax]/dt = k_InBax*Act*InBax - k_Bax*Bax - 1.0*Bax**4*k_o + 4*MAC*kr_o
    d[ActBcl2]/dt = k_BH3_Bcl2*Act*Bcl2 - kr_BH3Bcl2*ActBcl2
    d[EnaBcl2]/dt = k_BH3_Bcl2*Ena*Bcl2 - kr_BH3Bcl2*EnaBcl2
    d[MAC]/dt = 0.25*Bax**4*k_o - MAC*kr_o
    """

    p_name_map = {
        'one_step_BidT_BaxC_to_BidT_BaxA_kf': 'k_InBax',
        'reverse_BaxA_to_BaxC_k': 'k_Bax',
        'bind_BidT_Bcl2_kf': 'k_BH3_Bcl2',
        'bind_BidT_Bcl2_kr': 'kr_BH3Bcl2',
        'bind_BadM_Bcl2_kf': 'k_BH3_Bcl2',
        'bind_BadM_Bcl2_kr': 'kr_BH3Bcl2',
        'spontaneous_pore_BaxA_to_Bax4_kf': 'k_o',
        'spontaneous_pore_BaxA_to_Bax4_kr': 'kr_o'
    }

    s_name_map = {'s0': 'Act',
                  's1': 'Ena',
                  's2': 'InBax',
                  's3': 'Bcl2',
                  's4': 'Bax',
                  's5': 'ActBcl2',
                  's6': 'EnaBcl2',
                  's7': 'MAC'}

    generate_equations(model)

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
