from pysb import *
from local_macros import *
#from pysb.macros import catalyze_one_step_reversible
from shen_direct_module import direct_module

Model()

declare_monomers()

# MODEL TOPOLOGY
direct_module()

# Initial conditions
Initial(tBid(b=None), Parameter('tBid_0', 100))
Initial(Bax(b=None, state='C'), Parameter('Bax_0', 100))
Initial(Bcl2(b=None), Parameter('Bcl2_0', 100))
Initial(Bad(b=None), Parameter('Bad_0', 100))


