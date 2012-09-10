"""
Model M11a: Extrinsic apoptosis model incorporating the "Indirect" MOMP model
from [Chen2007febs]_.
"""

from pysb import *
from earm import shared
from earm.shared import V
from scipy.constants import N_A
from earm import albeck_modules
from earm import shen_modules
import re

Model()

# Declare monomers
albeck_modules.ligand_to_c8_monomers()
shen_modules.momp_monomers()
albeck_modules.apaf1_to_parp_monomers()

# The specific MOMP model to use
shen_modules.chen_febs_indirect(do_pore_assembly=True, do_pore_transport=True)

# Set initial condition for uncleaved Bid to 20nM, per the paper
Initial(Bid(state='U', bf=None), Parameter('Bid_0', 20e-9 * N_A * V))

albeck_modules.rec_to_bid()
albeck_modules.pore_to_parp()

# Declare common observables
shared.observables()

# Additional observables
Observable('aBax_', Bax(state='A', bf=None))
Observable('Bcl2_', Bcl2(bf=None))
Observable('Bcl2_Bid_', Bcl2(bf=1) % Bid(bf=1))
Observable('Bcl2_Bax_', Bcl2(bf=1) % Bax(bf=1))

