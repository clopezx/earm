"""
TODO: Docstring
'embedded together' w earm1.0"""

from pysb import *
from earm2 import macros
from earm2 import earm2_modules
from earm2 import albeck_modules

Model()

earm2_modules.momp_monomers()
Observable('aBax', Bax(state='A'))
Observable('cSmac', Smac(state='A'))

# The specific MOMP model to use
earm2_modules.direct()

# Set Bid initial condition to be tBid
model.update_initial_condition_pattern(Bid(state='U', bf=None),
                                       Bid(state='T', bf=None))
