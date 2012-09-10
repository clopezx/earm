"""
Model M8b: MOMP model "Current model + cooperativity" (Figure 11f) from [Albeck2008]_.
"""

from pysb import *
from earm import albeck_modules

Model()

albeck_modules.momp_monomers()

# The specific MOMP model to use
albeck_modules.albeck_11f()

# Observables
Observable('AcBax_', Bax(bf=None, state='A'))
Observable('cSmac', Smac(state='C'))

