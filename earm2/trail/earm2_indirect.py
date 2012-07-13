"""'direct' w earm1.0"""

from pysb import *
from earm2 import macros
from earm2 import earm2_modules
from earm2 import albeck_modules

Model()

# Declare monomers
albeck_modules.ligand_to_c8_monomers()
earm2_modules.momp_monomers()
albeck_modules.apaf1_to_parp_monomers()

# Declare observables and initial conditions
macros.all_observables()
macros.all_initial_conditions()

# Generate the upstream and downstream sections from the EARM 1.0 model
albeck_modules.rec_to_bid()
albeck_modules.pore_to_parp()

# The specific MOMP model to use
earm2_modules.indirect()

