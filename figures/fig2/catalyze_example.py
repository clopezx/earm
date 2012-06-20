from pysb import *
from pysb.macros import catalyze
from pysb.bng import generate_network, generate_equations

Model()

Monomer('C8', ['bf'])
Monomer('Bid', ['bf', 'state'], {'state': ['C', 'T']})

klist = [Parameter('kf', 1), Parameter('kr', 1), Parameter('kc', 1)]

catalyze(C8, 'bf', Bid(state='C'), 'bf', Bid(state='T'), klist)

Initial(C8(bf=None), Parameter('C8_0', 100))
Initial(Bid(bf=None, state='C'), Parameter('Bid_0', 100))

print (generate_network(model))
generate_equations(model)
print (model.odes)
