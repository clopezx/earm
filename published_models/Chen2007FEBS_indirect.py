from pysb import *
from local_macros import *
#from pysb.macros import catalyze_one_step_reversible

Model()

declare_monomers()

# MODEL TOPOLOGY
bind_table( [[  tBid, Bax(state='C')],
     [Bcl2,  (1, 1),         (1, 1)]])

pore_assembly(Bax(state='A'), Pore(), 4, [1, 1])

# Initial conditions
Initial(tBid(b=None), Parameter('tBid_0', 100))
Initial(Bax(b=None, state='C'), Parameter('Bax_0', 100))
Initial(Bcl2(b=None), Parameter('Bcl2_0', 100))

