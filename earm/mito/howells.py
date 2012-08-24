"""
Model M15b: MOMP model from Howells et al. (2010) J. Theor. Biol.

Howells, C. C., Baumann, W. T., Samuels, D. C., & Finkielstein, C. V. (2010).
The Bcl-2-associated death promoter (BAD) lowers the threshold at which the
Bcl-2-interacting domain death agonist (BID) triggers mitochondria
disintegration. Journal of Theoretical Biology.
:doi:`10.1016/j.jtbi.2010.11.040` :pmid:`21130780`.
"""

from pysb import *
from earm import shen_modules

Model()

shen_modules.momp_monomers()

# The specific MOMP model to use
shen_modules.howells()

# Observables
Observable('AcBax_', Bax(bf=None, state='A'))
Observable('Bax4_', MatchOnce(Bax(bf=None, s1=1, s2=4) %
                              Bax(bf=None, s1=2, s2=1) %
                              Bax(bf=None, s1=3, s2=2) %
                              Bax(bf=None, s1=4, s2=3)))
Observable('Bid_', Bid(bf=None))
#Observable('Bad_', Bad(bf=None))
Observable('Bcl2_', Bcl2(bf=None))
Observable('Bad_Bcl2_', Bcl2(bf=1) % Bad(bf=1))
Observable('Bid_Bcl2_', Bcl2(bf=1) % Bid(bf=1))
Observable('Bax_Bcl2_', Bcl2(bf=1) % Bax(bf=1))
Observable('pBad1433_', Bad(bf=None, state='C', serine='B'))
