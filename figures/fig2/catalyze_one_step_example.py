from pysb import *
from pysb.macros import catalyze_one_step
from pysb.bng import generate_network, generate_equations

Model()

Monomer('C8', ['bf'])
Monomer('Bid', ['bf', 'state'], {'state': ['U', 'T']})

kf = Parameter('kf', 1)

catalyze_one_step(C8(bf=None), Bid(state='U', bf=None),
                  Bid(state='T', bf=None), kf)

Initial(C8(bf=None), Parameter('C8_0', 100))
Initial(Bid(bf=None, state='U'), Parameter('Bid_0', 100))

print (generate_network(model))
generate_equations(model)

print (zip(model.species, model.odes))
