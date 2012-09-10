"""
Model M5b: MOMP model "Model B + Bax multimerization" (Figure 11c) from
[Albeck2008]_.
"""

from pysb import *
from earm import albeck_modules

Model()

albeck_modules.momp_monomers()

# The specific MOMP model to use
albeck_modules.albeck_11c()

# Observables
Observable('AcBax_', Bax(bf=None, state='A'))
Observable('cSmac', Smac(state='C'))

