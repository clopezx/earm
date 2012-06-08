from pysb import *
from local_macros import *

"""
Model drawn from:

Chun Chen, Jun Cui, Haizhu Lu, Rui Wang, Shuai Zhang, and Pingping Shen (2007)
"Modeling of the role of a Bax-activation switch in the mitochondrial apoptosis
decision," Biophys J 92:4304-4315
"""

Model()

declare_monomers()

# MODEL TOPOLOGY
catalyze_one_step_reversible(tBid, Bax(state='C'), Bax(state='A'), [1, 1])

bind_table([[  tBid,  Bax(state='A')],
     [Bcl2,  (1, 1),          (1, 1)]])

pore_assembly(Bax(state='A'), Pore(), 4, [1, 1])

displace_reversibly(Bcl2, Bax(state='A'), tBid, [1, 1])


# Initial conditions
Initial(tBid(b=None), Parameter('tBid_0', 100))
Initial(Bax(b=None, state='C'), Parameter('Bax_0', 100))
Initial(Bcl2(b=None), Parameter('Bcl2_0', 100))

