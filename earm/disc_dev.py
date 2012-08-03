"""Novel mechanism for DISC along with lopez_embedded and albeck apaf1_to_parp modules"""

from pysb import *
from earm import macros
from earm import disc_modules_dev
from earm import lopez_modules
from earm import albeck_modules

Model()

# Declare monomers
disc_modules_dev.lig_to_bid_monomers()
lopez_modules.momp_monomers()
albeck_modules.apaf1_to_parp_monomers()

# instantiate the momp and parp modules to use
lopez_modules.embedded()
albeck_modules.pore_to_parp()

# load the initial concentrations for the DISC module
disc_modules_dev.declare_initial_conditions()

# load the rates for the DISC module
disc_modules_dev.lig_to_bid_rates()

# specific DISC model to use
disc_modules_dev.lig_to_bid()

# Declare observables
macros.shared_observables()

