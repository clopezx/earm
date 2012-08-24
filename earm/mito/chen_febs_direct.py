"""
Model M10b: "Direct" MOMP model from Chen et al. (2007) FEBS Letters.

Chen, C., Cui, J., Zhang, W., & Shen, P. (2007). Robustness analysis identifies
the plausible model of the Bcl-2 apoptotic switch. FEBS letters, 581(26),
5143-5150. :doi:`10.1016/j.febslet.2007.09.063` :pmid:`17936275`.
"""

from pysb import *
from earm import shen_modules

Model()

shen_modules.momp_monomers()

# The specific MOMP model to use
shen_modules.chen_febs_direct()

# Observables
Observable('AcBax_', Bax(bf=None, state='A'))
Observable('Bax4_', MatchOnce(Bax(bf=None, s1=1, s2=4) %
                              Bax(bf=None, s1=2, s2=1) %
                              Bax(bf=None, s1=3, s2=2) %
                              Bax(bf=None, s1=4, s2=3)))
Observable('Bid_', Bid(bf=None))
Observable('Bad_', Bad(bf=None))
Observable('Bcl2_', Bcl2(bf=None))
Observable('Bcl2_Bad_', Bcl2(bf=1) % Bad(bf=1))
Observable('Bcl2_Bid_', Bcl2(bf=1) % Bid(bf=1))
Observable('Bcl2_Bax_', Bcl2(bf=1) % Bax(bf=1))
