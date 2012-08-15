"""
A variant version of Lopez embedded with no Smac/Cyto C to simplify
visualizations.
'embedded together' w earm1.0"""

from pysb import *
from earm import macros
from earm import lopez_modules
from earm import albeck_modules

Model()

lopez_modules.momp_monomers()
Observable('aBax', Bax(state='A'))
Observable('cSmac', Smac(state='C'))

# The specific MOMP model to use
lopez_modules.embedded(do_pore_transport=False)

# Set Bid initial condition to be tBid
model.update_initial_condition_pattern(Bid(state='U', bf=None),
                                       Bid(state='T', bf=None))
# Get rid of CytoC and Smac
model.parameters['Smac_0'].value = 0
model.parameters['CytoC_0'].value = 0

# Put some Noxa at the membrane
Initial(NOXA(state='M', bf=None), Parameter('mNOXA_0', 1))
Initial(Bad(state='M', bf=None), Parameter('mBad_0', 1))
