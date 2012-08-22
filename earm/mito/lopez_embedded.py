"""
Model M1b: "Embedded together" model of MOMP.
"""

from pysb import *
from earm import lopez_modules
from earm import albeck_modules

Model()

lopez_modules.momp_monomers()
Observable('aBax', Bax(state='A'))
Observable('cSmac', Smac(state='C'))

# The specific MOMP model to use
lopez_modules.embedded()

# Set Bid initial condition to be tBid
model.update_initial_condition_pattern(Bid(state='U', bf=None),
                                       Bid(state='T', bf=None))


