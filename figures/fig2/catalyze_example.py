from pysb import *
#from pysb.macros import catalyze
from pysb.bng import generate_network, generate_equations


def catalyze(enz, e_site, sub, s_site, prod, klist):
    kf, kr, kc = klist   # Get the parameters from the list

    # Create the rules
    rb = Rule('bind_%s_%s' % (enz.monomer.name, sub.monomer.name),
           enz({e_site:None}) + sub({s_site:None}) <>
           enz({e_site:1}) % sub({s_site:1}),
           kf, kr)
    rc = Rule('catalyze_%s%s_to_%s' %
           (enz.monomer.name, sub.monomer.name, prod.monomer.name),
           enz({e_site:1}) % sub({s_site:1}) >>
           enz({e_site:None}) + prod({s_site:None}),
           kc)
    return [rb, rc]

Model()

Monomer('C8', ['bf'])
Monomer('Bid', ['bf', 'state'], {'state': ['C', 'T']})

klist = [Parameter('kf', 1), Parameter('kr', 1), Parameter('kc', 1)]

catalyze(C8(), 'bf', Bid(state='C'), 'bf', Bid(state='T'), klist)

Initial(C8(bf=None), Parameter('C8_0', 100))
Initial(Bid(bf=None, state='C'), Parameter('Bid_0', 100))

print (generate_network(model))
generate_equations(model)

print (zip(model.species, model.odes))
