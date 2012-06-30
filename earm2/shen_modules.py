from pysb import *
from macros import *
import pysb.macros as macros
from pysb.util import alias_model_components
#from pysb.macros import catalyze_one_step_reversible

# Nomenclature from the paper:
# Ena = Bad
# Act = tBid
# Bcl2 = Bcl2
# InBax = Bax(state='C')
# AcBax = Bax(state='A')
# MAC = Pore()

transloc_rates = [1e-2, 1e-2]

def chen2007BiophysJ():
    alias_model_components()

    # One-step activation of Bax by tBid
    macros.catalyze_one_step_reversible(
        Bid(state='T', bf=None), Bax(state='C', bf=None),
        Bax(state='A', bf=None), [1, 1])

    active_monomer = {'state': 'A', 's1': None, 's2': None}
    bind_table([[                        Bcl2],
                [Bid(state='T'),       (1, 1)],
                [Bax(active_monomer),  (1, 1)]])

    displace_reversibly(Bcl2, Bax(active_monomer), Bid(state='T'), [1, 1])

    #assemble_pore(Bax(state='A'), 4, Pore(),
    #              [Parameter('p1',1), Parameter('pr', 1)])

    #catalyze(Pore(), CytoC(state='M'), CytoC(state='C'), [1, 1, 1])

def howells2011():
    chen2007BiophysJ()
    # Reactions regarding Bad
    #bind(Bcl2, Bad, [1, 1])
    #displace(tBid(s=1) % Bcl2(s=1), Bad, Bad(s=1) % Bcl2(s = 1), k)
    #two_state_equilibrium(Bad(loc='M'), Bad(loc='C'), [kf, kr])
    #transition_table([[Bad(s=1) % Bcl2(s=1), Bad()????
    #                  [pBad, pBad:14-3-3, k]
    #                  [pBad:14-3-3, Bad, kBadRel],
    #                  [Bad, pBad, kBADphos1]])

def chen2007FEBS_indirect():
    alias_model_components()
    bind_table([[         tBid,  Bax(state='C')],
                [Bcl2,  (1, 1),          (1, 1)]])
    assemble_pore(Bax(state='C'), Pore(), 4, [1, 1])

def chen2007FEBS_direct():
    alias_model_components()
    catalyze_one_step_reversible(tBid, Bax(state='C'), Bax(state='A'),
                                [1, 1])
    bind_table([[         tBid,     Bad],
                [Bcl2,  (1, 1),  (1, 1)]])
    assemble_pore(Bax(state='A'), Pore(), 4, [1, 1])

def cui2008_direct():
    alias_model_components()
    chen2007FEBS_direct()
    synthesize_degrade_table(
       [[Bax(state='C'),                   0.06,  0.001],
        [Bax(state='A'),                   None,  0.001],                     
        [tBid,                            0.001,  0.001],
        [Bcl2,                             0.03,  0.001],
        [tBid(b=1) % Bcl2(b=1),            None,  0.005],
        [Bax(b=1, state='C') % Bcl2(b=1),  None,  0.005],
        [Bad,                             0.001,  0.001],
        [Bad(b=1) % Bcl2(b=1),             None,  0.005],
        [Pore,                             None, 0.0005]])

def cui2008_direct1():
    alias_model_components()
    cui2008_direct()
    bind(Bax(state='A'), Bcl2, [1, 1])

def cui2008_direct2():
    alias_model_components()
    cui2008_direct1()
    catalyze(Bax(state='A'), Bax(state='C'), Bax(state='A'), [1, 1, 1])




