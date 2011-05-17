from pysb import *
from pysbhelperfuncs import *
# Show how to build a complex model
# Based on Letai data

Model()
import letai_parms # Parameters

# Species monomer sections
Monomer('tBid',    ['bf', 'state'], {'state':['A', 'I']})
Monomer('Bax',     ['bf', 'state'], {'state':['A', 'I']})
Monomer('Bak',     ['bf', 'state'], {'state':['A', 'I']})
Monomer('CytC',    ['bf', 'state'], {'state':['cyto', 'mito']})
Monomer('BaxPore', ['bf'])
Monomer('BakPore', ['bf'])
Monomer('Bcl2',    ['bf'])
Monomer('Bad',     ['bf'])
# Initial species
Initial(tBid(bf = None, state = 'A'), tBid_0)
Initial(Bax(bf=None, state='I'), Bax_0)
Initial(CytC(bf=None, state='mito'), CytC_0)
Initial(Bcl2(bf=None), Bcl2_0)
Initial(Bad(bf=None), Bad_0)

# Rules

# Bax binds to active tBid, Bax becomes Active
twostepact(tBid(state='A'), Bax(state='I'), Bax(bf = None, state='A'),
           kbidbaxf, kbidbaxr, kbidbaxc)

# Active Bax and Bak aggregates into a pore
twostepconv(Bax(state='A'), Bax(state='A'), BaxPore(bf=None),
           kbaxporef, kbaxporer, kbaxporec)
twostepconv(Bak(state='A'), Bak(state='A'), BakPore(bf=None),
           kbakporef, kbakporer, kbakporec)

# Active Bax and Bak pores transport CytC out of the membrane in two-step rxns
twostepact(BaxPore(), CytC(state='cyto'), CytC(bf=None, state='mito'),
           kbaxcytcf, kbaxcytcr, kbaxcytcc)
twostepact(BakPore(), CytC(state='cyto'), CytC(bf=None, state='mito'),
           kbakcytcf, kbakcytcr, kbakcytcc)

# Bcl2 inhibitors of Bax and Bak
Inhtable([[        Bcl2, BclxL, Mcl1, Bclw, Bfl1]
          [Bax, baxbcl2, baxbclxl, baxmcl1, ]
          [Bak                              ]])



# Bcl-2 inhibits Bax
simplebind(Bax(state='A'), Bcl2(), kbaxbcl2f, kbaxbcl2r)

# Bad inhibits Bcl-2
simplebind(Bcl2(), Bad(), kbcl2badf, kbcl2badr)
