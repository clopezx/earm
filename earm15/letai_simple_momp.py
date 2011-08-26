from pysb import *
from pysbhelperfuncs import *
# Show how to incorporate complex rules easily
# Based on Letai data

Model()

# Species monomer sections
Monomer('tBid',    ['bf', 'state'], {'state':['A', 'I']})
Monomer('Bax',     ['bf', 'state'], {'state':['A', 'I']})
Monomer('CytC',    ['bf', 'state'], {'state':['cyto', 'mito']})
Monomer('BaxPore', ['bf'])
Monomer('Bcl2',    ['bf'])
Monomer('Bad',     ['bf'])

# Parameter section
# Reaction rates
Parameter('kbidbaxf',  9.99)
Parameter('kbidbaxr',  9.99)
Parameter('kbidbaxc',  9.99)
Parameter('kbaxporef', 8.88)
Parameter('kbaxporer', 8.88)
Parameter('kbaxporec', 8.88)
Parameter('kbaxcytcf', 7.77)
Parameter('kbaxcytcr', 7.77)
Parameter('kbaxcytcc', 7,77)
Parameter('kbaxbcl2f', 6.66)
Parameter('kbaxbcl2r', 6.66)
Parameter('kbcl2badf', 5.55)
Parameter('kbcl2badr', 5.55)
# Initial conditions
Parameter('tBid_0', 1000)
Parameter('Bax_0',  1000) 
Parameter('CytC_0', 1000) 
Parameter('BaxPore_0', 0.0)
Parameter('Bcl2_0', 1000) 
Parameter('Bad_0',  1000) 
# Initial species
Initial(tBid(bf = None, state = 'A'), tBid_0)
Initial(Bax(bf=None, state='I'), Bax_0)
Initial(CytC(bf=None, state='mito'), CytC_0)
Initial(BaxPore(bf=None), BaxPore_0)
Initial(Bcl2(bf=None), Bcl2_0)
Initial(Bad(bf=None), Bad_0)

# Rules

# Bax binds to active tBid, Bax becomes Active
twostepact(tBid(state='A'), Bax(state='I'), Bax(bf = None, state='A'),
           kbidbaxf, kbidbaxr, kbidbaxc)

# Active Bax aggregates into a pore
twostepconv(Bax(state='A'), Bax(state='A'), BaxPore(bf=None),
           kbaxporef, kbaxporer, kbaxporec)

# Active pore can transport CytC out of the membrane in two-step rxn
twostepact(BaxPore(), CytC(state='cyto'), CytC(bf=None, state='mito'),
           kbaxcytcf, kbaxcytcr, kbaxcytcc)

# Bcl-2 inhibits Bax
simplebind(Bax(state='A'), Bcl2(), kbaxbcl2f, kbaxbcl2r)

# Bad inhibits Bcl-2
simplebind(Bcl2(), Bad(), kbcl2badf, kbcl2badr)
