from pysb import *
from pysbhelperfuncs import *
# Show how to incorporate complex rules easily
# Based on Letai data

Model()

# Species monomer sections
Monomer('tBid', ['bf', 'state'], {'state':['A', 'I']})
Monomer('Bax', ['bf', 'state'], {'state':['A', 'I']})
Monomer('CytC', ['b', 'state'], {'state':['mito', 'cyto']})
Monomer('BaxPore', ['b'])

# Parameter section
Parameter('kbidbaxf',  9.99)
Parameter('kbidbaxr',  9.99)
Parameter('kbidbaxc',  9.99)
Parameter('kbaxporef', 8.88)
Parameter('kbaxporer', 8.88)
Parameter('kbaxporec', 8.88)

# Bax binds to active tBid, Bax becomes Active
catact(tBid(state='A'), Bax(state='I'), Bax(bf = None, state='A'),
       kbidbaxf, kbidbaxr, kbidbaxc)

# Active Bax aggregates into a pore
catact(Bax(state='A'), Bax(state='A'), BaxPore(b=None),
       kbaxporef, kbaxporer, kbaxporec)
