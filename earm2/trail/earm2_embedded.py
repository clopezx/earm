"""'embedded together' w earm1.0"""

from pysb import *
from . import macros
from . import earm2_modules
from . import albeck_modules

Model()

macros.declare_all_monomers()

earm2_modules.embedded()

# Generate the upstream and downstream sections from the EARM 1.0 model
albeck_modules.rec_to_bid()
albeck_modules.pore_to_parp()

macros.declare_all_initial_conditions()

# Observables
# TODO: figure out where these should live
Observe('mBid',  Bid(state='M'))
Observe('cSmac', Smac(state='A'))
Observe('cPARP', PARP(state='C'))
