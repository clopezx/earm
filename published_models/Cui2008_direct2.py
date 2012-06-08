from pysb import *
from local_macros import *
from shen_direct_module import direct_module

from Cui2008_direct import model as direct_model

Model(base=direct_model)

# MODEL TOPOLOGY
bind(Bax(state='A'), Bcl2, [1, 1])
catalyze(Bax(state='A'), Bax(state='C'), Bax(state='A'), [1, 1, 1])

