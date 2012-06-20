from pysb import *
from pysb.macros import pore_assembly
from pysb.bng import generate_network, generate_equations

Model()

Monomer('Bax', ['bf', 'bh3', 'd2'])

krow = [Parameter('kf', 1), Parameter('kr', 1)]
klist = [krow for i in range(4)]

pore_assembly(Bax, 4, klist)

Initial(Bax(bf=None, s1=None, s2=None), Parameter('Bax_0', 100))

print (generate_network(model))
generate_equations(model)
print (model.odes)


