from pysb import *
from local_macros import *

Model()

catalyze(tBid, Bax(state='C'), Bax(state='A'), [1, 1, 1])

