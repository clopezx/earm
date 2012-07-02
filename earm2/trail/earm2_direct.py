"""'direct' w earm1.0"""

from pysb import *
from earm2 import macros
from earm2 import earm2_modules
from earm2 import albeck_modules

Model()

macros.all_monomers()
macros.all_observables()

# Generate the upstream and downstream sections from the EARM 1.0 model
albeck_modules.rec_to_bid()
albeck_modules.pore_to_parp()

# The specific MOMP model to use
earm2_modules.direct()
macros.all_initial_conditions()

