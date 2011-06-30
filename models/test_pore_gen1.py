#!/usr/bin/env python2.6

from pysb import *
from pysbhelperfuncs import *

Model()

pore_size = 5

pore_params = []
for k in range(1, pore_size+1):
    pore_params.append([
            Parameter('kpore%df' % k, k*10),
            Parameter('kpore%dr' % k, k)
            ])
Parameter('Bak_0', 1e3)
    
Monomer('Bak', 'b')

pore_assembly(Bak(), 'b', pore_size, pore_params)

Initial(Bak(b=None), Bak_0)
