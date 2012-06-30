"""Model from Chen 2007, Biophys J."""

from pysb import *
from earm2 import macros
from earm2 import earm2_modules
from earm2 import albeck_modules
from earm2 import shen_modules

Model()

macros.momp_monomers()
Observable('aBax', Bax(state='A'))
Observable('cSmac', Smac(state='A'))

Initial(Bid(state='T', bf=None), Parameter('Bid_0', 100))
Initial(Bax(state='C', bf=None, s1=None, s2=None), Parameter('Bax_0', 100))
Initial(Bcl2(bf=None), Parameter('Bcl2_0', 100))

# The specific MOMP model to use
shen_modules.chen2007BiophysJ()



