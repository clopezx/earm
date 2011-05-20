from pysb import *
from pysbhelperfuncs import *
# Show how to build a complex model
# Based on Letai data

Model()
#FIXME! import letai_parms # Parameters

# Monomers 
Monomer('tBid',    ['bf', 'state'], {'state':['A', 'I']})
Monomer('Bax',     ['bf', 'state'], {'state':['A', 'I']})
Monomer('Bak',     ['bf', 'state'], {'state':['A', 'I']})
Monomer('CytC',    ['bf', 'state'], {'state':['cyto', 'mito']})
Monomer('BaxPore', ['bf'])
Monomer('BakPore', ['bf'])
Monomer('Bcl2',    ['bf'])
Monomer('BclxL',   ['bf'])
Monomer('Mcl1',    ['bf'])
Monomer('Bad',     ['bf'])

# Reaction rates
Parameter('kbidbaxf',  9.99)
Parameter('kbidbaxr',  9.99)
Parameter('kbidbaxc',  9.99)
Parameter('kbaxporef', 8.88)
Parameter('kbaxporer', 8.88)
Parameter('kbaxporec', 8.88)
Parameter('kbakporef', 8.85)
Parameter('kbakporer', 8.85)
Parameter('kbakporec', 8.85)
Parameter('kbaxcytcf', 7.77)
Parameter('kbaxcytcr', 7.77)
Parameter('kbaxcytcc', 7,77)
Parameter('kbakcytcf', 7.75)
Parameter('kbakcytcr', 7.75)
Parameter('kbakcytcc', 7,75)
Parameter('kbaxbcl2f', 6.66)
Parameter('kbaxbcl2r', 6.66)
Parameter('kbcl2badf', 5.55)
Parameter('kbcl2badr', 5.55)
#---inhibitors parameters---
Parameter('baxbcl2f', 3.33)
Parameter('baxbcl2r', 3.33)
Parameter('baxbclxlf', 3.33)
Parameter('baxbclxlr', 3.33)
Parameter('baxmcl1f', 3.33)
Parameter('baxmcl1r', 3.33)
Parameter('bakbcl2f', 3.33)
Parameter('bakbcl2r', 3.33)
Parameter('bakbclxlf', 3.33)
Parameter('bakbclxlr', 3.33)
Parameter('bakmcl1f', 3.33)
Parameter('bakmcl1r', 3.33)


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
sbindtable([[                     Bcl2, BclxL,  Mcl1],
            [                       {},    {},    {}],
            [Bax, {'state':'A'},  True,  True, False],
            [Bak, {'state':'A'}, False,  True,  True]],
           model.parameters)

# Bad inhibits Bcl-2
simplebind(Bcl2(), Bad(), kbcl2badf, kbcl2badr)


# Initial conditions
# Initial concentrations
Parameter('tBid_0', 1000)
Parameter('Bax_0',  1000) 
Parameter('Bak_0',  1000) 
Parameter('CytC_0', 1000) 
Parameter('BaxPore_0', 0.0)
Parameter('BakPore_0', 0.0)
Parameter('Bcl2_0', 1000) 
Parameter('BclxL_0', 1000) 
Parameter('Mcl1_0',  1000) 
Parameter('Bad_0', 1000)
# Initial species
Initial(tBid(bf = None, state = 'A'), tBid_0)
Initial(Bax(bf=None, state='I'), Bax_0)
Initial(CytC(bf=None, state='mito'), CytC_0)
Initial(Bcl2(bf=None), Bcl2_0)
Initial(Bad(bf=None), Bad_0)
Initial(Bak(bf =None, state='I'), Bak_0)
Initial(BaxPore(bf=None), BaxPore_0)
Initial(BakPore(bf=None), BakPore_0)
Initial(BclxL(bf=None), BclxL_0)
Initial(Mcl1(bf=None), Mcl1_0)

