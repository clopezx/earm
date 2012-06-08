from pysb import *
from local_macros import *

Model()

declare_monomers()

catalyze(tBid, Bax(state='C'), Bax(state='A'), [1, 1, 1])

bind_table([ [Bcl2,    Mcl1],
     [tBid,  (1,1),   (1,1)]])

displace_reversibly(Bcl2, tBid, Bad, [1, 1])

displace_reversibly(Mcl1, tBid, Noxa, [1, 1])


