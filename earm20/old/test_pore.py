#!/usr/bin/env python2.6

from pysb import *
from pysbhelperfuncs import *

Model()

v = .07; # mitochondria compartment volume/cell volume
max_pore_size = 5
min_pore_transport_size = 4

Monomer('Bak', ['bh3', 'd2', 'bf'])
Monomer('Mcl1', ['bh3'])
Monomer('Smac', ['bf', 'loc'], {'loc': ['M','C']})

Parameter('Bak_0', 1e5)
Parameter('Mcl1_0', 2e4)
Parameter('Smac_0', 1e5)
    
Initial(Bak(bh3=None, d2=None, bf=None), Bak_0)
Initial(Mcl1(bh3=None), Mcl1_0)
Initial(Smac(bf=None, loc='M'), Smac_0)

pore_assembly_params = []
for k in range(1, max_pore_size+1):
    pore_assembly_params.append([
            Parameter('kbakpore%df' % k, 1e-6/v*2),
            Parameter('kbakpore%dr' % k, 1e-3)
            ])
pore_assembly(Bak(bf=None), max_pore_size, pore_assembly_params)

pore_transport_params = []
for k in range(min_pore_transport_size, max_pore_size+1):
    pore_transport_params.append([
            Parameter('ksmactrans%df' % k, 2e-6/v),
            Parameter('ksmactrans%dr' % k, 1e-3),
            Parameter('ksmactrans%dc' % k, 1e+1)
            ])
pore_transport(Bak(bf=None), Smac(loc='M'), Smac(loc='C'),
               min_pore_transport_size, max_pore_size, pore_transport_params)

Parameter('kmcl1bakf', 1.0)
Parameter('kmcl1bakr', 3.0e-1)
Rule('mcl1_inhibit_bak', Bak(bh3=None, d2=None) + Mcl1(bh3=None) <> Bak(bh3=1, d2=None) % Mcl1(bh3=1),
     kmcl1bakf, kmcl1bakr)
