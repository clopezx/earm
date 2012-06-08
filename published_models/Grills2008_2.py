from pysb import *
from local_macros import *

Model()

declare_monomers()

catalyze(tBid, Bax(state='C'), Bax(state='A'), [1, 1, 1])

displace_reversibly(Bcl2, tBid, Bad, [1, 1])

