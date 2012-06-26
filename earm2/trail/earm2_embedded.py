"""'embedded together' w earm1.0"""

from pysb import *
from earm2 import macros
from earm2 import earm2_modules
from earm2 import albeck_modules

Model()

macros.declare_all_monomers()

earm2_modules.embedded()

# Generate the upstream and downstream sections from the EARM 1.0 model
albeck_modules.rec_to_bid()
albeck_modules.pore_to_parp()

macros.declare_all_initial_conditions()

# Observables
# TODO: figure out where these should live
Observable('mBid',  Bid(state='M'))
Observable('cSmac', Smac(state='A'))
Observable('cPARP', PARP(state='C'))
