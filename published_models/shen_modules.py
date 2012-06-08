from pysb import *
from local_macros import *
import pysb.macros as macros
#from pysb.macros import catalyze_one_step_reversible

# Nomenclature from the paper:
# Ena = Bad
# Act = tBid
# Bcl2 = Bcl2
# InBax = Bax(state='C')
# AcBax = Bax(state='A')
# MAC = Pore()


def chen2007BiophysJ():
    macros.alias_model_components()

    catalyze_one_step_reversible(tBid, Bax(state='C'), Bax(state='A'), [1, 1])

    bind_table([[         tBid,  Bax(state='A')],
                [Bcl2,  (1, 1),         (1, 1)]])

    assemble_pore(Bax(state='A'), Pore(), 4, [1, 1])

    displace_reversibly(Bcl2, Bax(state='A'), tBid, [1, 1])


def chen2007FEBS_indirect():
    macros.alias_model_components()

    bind_table([[         tBid,  Bax(state='C')],
                [Bcl2,  (1, 1),          (1, 1)]])

    assemble_pore(Bax(state='A'), Pore(), 4, [1, 1])


def chen2007FEBS_direct():
    macros.alias_model_components()

    catalyze_one_step_reversible(tBid, Bax(state='C'), Bax(state='A'),
                                [1, 1])

    bind_table([[         tBid,     Bad],
                [Bcl2,  (1, 1),  (1, 1)]])

    assemble_pore(Bax(state='A'), Pore(), 4, [1, 1])


def cui2008_direct():
    macros.alias_model_components()

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
    macros.alias_model_components()

    cui2008_direct()

    bind(Bax(state='A'), Bcl2, [1, 1])


def cui2008_direct2():
    macros.alias_model_components()

    cui2008_direct1()

    catalyze(Bax(state='A'), Bax(state='C'), Bax(state='A'), [1, 1, 1])


