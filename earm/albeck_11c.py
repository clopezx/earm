"""
Model M5a: Extrinsic apoptosis model incorporating the MOMP model "Model B +
Bax multimerization" (Figure 11c) from Albeck et al. (2008) PLoS Biology.

Albeck, J. G., Burke, J. M., Spencer, S. L., Lauffenburger, D. A., and
Sorger, P. K. (2008). Modeling a snap-action, variable-delay switch
controlling extrinsic cell death. PLoS Biology, 6(12), 2831-2852.
:doi:`10.1371/journal.pbio.0060299` :pmid:`19053173`.
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
albeck_modules.albeck_11c()

# Declare shared observables
shared.observables()

# Additional observables
Observable('aBax_', Bax(bf=None, state='A'))

