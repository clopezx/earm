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
shen_modules.chen2007FEBS_indirect(pore_assembly=True)

# Observables
Observable('Bax4_', MatchOnce(Bax(bf=None, s1=1, s2=4) %
                              Bax(bf=None, s1=2, s2=1) %
                              Bax(bf=None, s1=3, s2=2) %
                              Bax(bf=None, s1=4, s2=3)))
Observable('Bid_', Bid(bf=None))
Observable('Bcl2_', Bcl2(bf=None))
Observable('Bcl2_Bid_', Bcl2(bf=1) % Bid(bf=1))
Observable('Bcl2_Bax_', Bcl2(bf=1) % Bax(bf=1))


def print_original_odes():
    """Show the ODEs using the nomenclature from the paper.

        >>> import earm2.mito.chen2007FEBS_indirect as i
        >>> i.print_original_odes()
        d[BH3]/dt = -k_BH3_Bcl2*BH3*Bcl2 + kr_BH3Bcl2*BH3Bcl2
        d[Bax]/dt = -k_Bax_Bcl2*Bax*Bcl2 + kr_BaxBcl2*BaxBcl2 - 1.0*Bax**4*k_o + 4*MAC*kr_o
        d[Bcl2]/dt = -k_Bax_Bcl2*Bax*Bcl2 + kr_BaxBcl2*BaxBcl2 - k_BH3_Bcl2*BH3*Bcl2 + kr_BH3Bcl2*BH3Bcl2
        d[BH3Bcl2]/dt = k_BH3_Bcl2*BH3*Bcl2 - kr_BH3Bcl2*BH3Bcl2
        d[BaxBcl2]/dt = k_Bax_Bcl2*Bax*Bcl2 - kr_BaxBcl2*BaxBcl2
        d[MAC]/dt = 0.25*Bax**4*k_o - MAC*kr_o

    """

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

if __name__ == "__main__":
    import doctest
    doctest.testmod()
