"""
Model M7b: MOMP model "Current model" (Figure 11e) from [Albeck2008]_.
"""

from pysb import *
from earm import albeck_modules

Model()

albeck_modules.momp_monomers()

# The specific MOMP model to use
albeck_modules.albeck_11e()

# Observables
Observable('AcBax_', Bax(bf=None, state='A'))
Observable('cSmac', Smac(state='C'))

