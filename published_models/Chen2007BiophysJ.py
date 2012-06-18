from pysb import *
from local_macros import *
from earm_2_0_macros import declare_monomers
import earm_1_0modules as earm_modules
import shen_modules

"""
Model drawn from:

Chun Chen, Jun Cui, Haizhu Lu, Rui Wang, Shuai Zhang, and Pingping Shen (2007)
"Modeling of the role of a Bax-activation switch in the mitochondrial apoptosis
decision," Biophys J 92:4304-4315
"""

Model()

declare_monomers()

from earm_2_0_parms import parameter_dict as kd

# Initial conditions
#Initial(tBid(b=None), Parameter('tBid_0', 100))
#Initial(Bax(b=None, state='C'), Parameter('Bax_0', 100))
#Initial(Bcl2(b=None), Parameter('Bcl2_0', 100))

earm_modules.rec_to_bid(kd)

shen_modules.chen2007BiophysJ()

earm_modules.pore_to_parp(kd)
