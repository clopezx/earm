"""'embedded together' w earm1.0"""

from pysb import *
from earm2 import macros
from earm2 import earm2_modules
from earm2 import albeck_modules

Model()

macros.momp_monomers()
Observable('aBax', Bax(state='A'))
Observable('cSmac', Smac(state='A'))

# The specific MOMP model to use
earm2_modules.embedded()
macros.momp_initial_conditions('embedded', bid_state='T')


