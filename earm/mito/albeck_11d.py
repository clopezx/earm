"""
Model M6b: MOMP model "Model C + mitochondrial transport" (Figure 11d) from
Albeck et al.  (2008) PLoS Biology.

Albeck, J. G., Burke, J. M., Spencer, S. L., Lauffenburger, D. A., and
Sorger, P. K. (2008). Modeling a snap-action, variable-delay switch
controlling extrinsic cell death. PLoS Biology, 6(12), 2831-2852.
:doi:`10.1371/journal.pbio.0060299` :pmid:`19053173`.
"""

from pysb import *
from earm import albeck_modules

Model()

albeck_modules.momp_monomers()

# The specific MOMP model to use
albeck_modules.albeck_11d()

# Observables
Observable('AcBax_', Bax(bf=None, state='A'))
Observable('cSmac', Smac(state='C'))

