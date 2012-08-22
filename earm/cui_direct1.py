"""
Model M13a: Extrinsic apoptosis model incorporating the "Direct 1" MOMP model
from Cui et al. (2008) PLoS One.

Cui, J., Chen, C., Lu, H., Sun, T., & Shen, P. (2008). Two independent positive
feedbacks and bistability in the Bcl-2 apoptotic switch. PLoS ONE, 3(1), e1469.
:doi:`10.1371/journal.pone.0001469` :pmid:`18213378`.
"""

from pysb import *
from earm import shared
from earm.shared import cell_vol
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
shen_modules.cui_direct1(do_pore_transport=True)

# Set initial condition for uncleaved Bid to 20nM, per the paper
Initial(Bid(state='U', bf=None), Parameter('Bid_0', 20e-9 * N_A * cell_vol))

albeck_modules.rec_to_bid()
albeck_modules.pore_to_parp()

# Declare common observables
shared.observables()

# Additional observables
Observable('aBax_', Bax(state='A', bf=None))
Observable('Bcl2_', Bcl2(bf=None))
Observable('Bcl2_Bid_', Bcl2(bf=1) % Bid(bf=1))
Observable('Bcl2_Bax_', Bcl2(bf=1) % Bax(bf=1))

