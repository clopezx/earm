from pysb import *
from earm2 import albeck_modules

Model()

albeck_modules.momp_monomers()

# The specific MOMP model to use
albeck_modules.albeck_11b(pore_transport=True)

# Observables
Observable('AcBax_', Bax(bf=None, state='A'))
Observable('Bid_', Bid(bf=None))
#Observable('Bcl2_', Bcl2(bf=None))
#Observable('Bcl2_Bid_', Bcl2(bf=1) % Bid(bf=1))
Observable('Bcl2_Bax_', Bcl2(bf=1) % Bax(bf=1))
Observable('cSmac', Smac(state='C'))

