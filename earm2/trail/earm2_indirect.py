"""'indirect' w earm1.0"""

from pysb import *
from earm2 import macros
from earm2 import earm2_modules
from earm2 import albeck_modules

Model()

macros.declare_all_monomers()
macros.declare_all_observables()

# Generate the upstream and downstream sections from the EARM 1.0 model
albeck_modules.rec_to_bid()
albeck_modules.pore_to_parp()

# The specific MOMP model to use
earm2_modules.indirect()
macros.declare_all_initial_conditions('indirect')


