from pysb import *
from local_macros import *
#from pysb.macros import catalyze_one_step_reversible

Model()

# Nomenclature from the paper:
# Ena = Bad
# Act = tBid
# Bcl2 = Bcl2
# InBax = Bax(state='C')
# AcBax = Bax(state='A')
# MAC = Pore()

declare_monomers()

catalyze_one_step_reversible(tBid, Bax(state='C'), Bax(state='A'), [1, 1, 1])

bind_table([[  tBid,     Bad],
     [Bcl2,  (1, 1),  (1, 1)]])

pore_assembly(Bax(state='A'), Pore(), 4, [1, 1])

