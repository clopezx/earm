from pysb import *
from local_macros import *
#from pysb.macros import catalyze_one_step_reversible

Model()

declare_monomers()

bind_table([[  tBid, Bax(state='C')],
     [Bcl2,  (1, 1),         (1, 1)]])

pore_assembly(Bax(state='A'), Pore(), 4, [1, 1])

