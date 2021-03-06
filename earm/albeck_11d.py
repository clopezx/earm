"""
Model M6a: Extrinsic apoptosis model incorporating the MOMP model "Model C +
mitochondrial transport" (Figure 11d) from [Albeck2008]_.
"""

from pysb import *
from earm import albeck_modules
from earm import shared

Model()

# Declare monomers
albeck_modules.all_monomers()

# Generate the upstream and downstream sections
albeck_modules.rec_to_bid()
albeck_modules.pore_to_parp()

# The specific MOMP model to use
albeck_modules.albeck_11d()

# Declare shared observables
shared.observables()

# Additional observables
Observable('aBax_', Bax(bf=None, state='A'))

