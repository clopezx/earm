from pysb.macros import *
from pysb import *
from pysb.bng import generate_network

Model()

Monomer('subunit', ['s1', 's2'])

assemble_pore_sequential(subunit, 's1', 's2', 4, [[1, 1], [1, 1], [1, 1]])

Initial(subunit(s1=None, s2=None), Parameter('subunit_0', 1))

Observable('sub', subunit())

generate_network(model)
