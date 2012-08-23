"""
Model M12b: "Direct" MOMP model from Cui et al. (2008) PLoS One.

Cui, J., Chen, C., Lu, H., Sun, T., & Shen, P. (2008). Two independent positive
feedbacks and bistability in the Bcl-2 apoptotic switch. PLoS ONE, 3(1), e1469.
:doi:`10.1371/journal.pone.0001469` :pmid:`18213378`.
"""

from pysb import *
from earm import shen_modules

Model()

shen_modules.momp_monomers()

# The specific MOMP model to use
shen_modules.cui_direct()

# Observables
Observable('AcBax_', Bax(bf=None, state='A'))
#Observable('Bax4_', MatchOnce(Bax(bf=None, s1=1, s2=4) %
#                              Bax(bf=None, s1=2, s2=1) %
#                              Bax(bf=None, s1=3, s2=2) %
#                              Bax(bf=None, s1=4, s2=3)))
Observable('Bid_', Bid(bf=None))
Observable('Bad_', Bad(bf=None))
Observable('Bcl2_', Bcl2(bf=None))
Observable('Bcl2_Bad_', Bcl2(bf=1) % Bad(bf=1))
Observable('Bcl2_Bid_', Bcl2(bf=1) % Bid(bf=1))
Observable('Bcl2_Bax_', Bcl2(bf=1) % Bax(bf=1))

