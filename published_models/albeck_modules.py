from pysb import *
from local_macros import *
import pysb.macros as macros

# From  PLoS biology paper

def albeck2008_11b():
    catalyze(tBid, Bax(state='C'), Bax(state='A'), [kf, kr, kc])
    bind(Bax(state='A'), Bcl2, [kf, kr])
    catalyze(Bax(state='A'), Smac(loc='c'), Smac(loc='r'), [kf, kr, kc])
    
def albeck2008_11c():
    catalyze(tBid, Bax(state='C'), Bax(state='A'), [kf, kr, kc])

    dimerize(Bax(state='A'), Bax2, [kf, kr])
    dimerize(Bax2, Bax4, [kf, kr])

    bind_table([[            Bax,      Bax2,     Bax4],
                [Bcl2,  [kf, kr],  [kf, kr],  [kf, kr]])

    catalyze(Bax4, Smac(loc='c'), Smac(loc='r'), [kf, kr, kc])

# Needs separate mitochondrial compartment--by rate scaling?
def albeck2008_11d():
    catalyze(tBid, Bax(state='C'), Bax(state='A'), [kf, kr, kc])

    dimerize(Bax(state='A'), Bax2, [kf, kr])
    dimerize(Bax2, Bax4, [kf, kr])

    bind_table([[            Bax,      Bax2,     Bax4],
                [Bcl2,  [kf, kr],  [kf, kr],  [kf, kr]])

    catalyze(Bax4, Smac(loc='c'), Smac(loc='r'), [kf, kr, kc])
        
def albeck2008_11e():
    catalyze(tBid, Bax(state='C'), Bax(state='A'), [kf, kr, kc])

    dimerize(Bax(state='A'), Bax2, [kf, kr])
    dimerize(Bax2, Bax4, [kf, kr])

    bind_and_convert(Bax4, M, Pore, [kf, kr])

    bind_table([[            Bax,      Bax2,     Bax4],
                [Bcl2,  (kf, kr),  (kf, kr),  (kf, kr)]])

    catalyze(Pore, Smac(loc='c'), Smac(loc='r'), [kf, kr, kc])


def albeck2008_11f():
    pass



