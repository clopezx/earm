from pysb import *

# Basic model aiming to build a model based on the data by 
# Andrews, Chen, and Kinnally

Model()

# Compartments section
Compartment('soln', dimension = 3, size = XX, parent = None)
Compartment('lipo_brane', dimension = 2, size = XX, parent = soln)
Compartment('liposome', dimension = 3, size = XX, parent = lipo_brane)

# Species monomer sections
Monomer('tBid', ['b1','b2', 'state'], {'state':['A', 'inA']})
Monomer('Bax', ['b', 'm1', 'm2', 'state'], {'state':['A', 'inA']})
Monomer('Cc-m', ['b'])

Parameter('[NAME]', [VALUE])

# tBid diffuses to the liposome membrane
Rule('tBid_to_brane',
     tBid(b1 = None, b2=None, s='A') ** soln <> tBid(b = None, b2=None, s='A') ** lipo_brane,
     ktbidbranef, ktbidbraner)

# Bax binds to a tBid on the membrane, Bax becomes Aive
Rule('Bax_to_Bidbrane',
     Bax(b = None, m1 = None, m2 = None, state='I') ** soln + 
     tBid(b1 = None, b2 = None, s='A') ** lipo_brane <>
     Bax(b = 1, m1 = None, m2 = None, state='I') ** lipo_brane % 
     tBid(b1 = 1, b2 = None s='A') ** lipo_brane,
     kbaxbranef, baxbraner)

# Active Bax brings another Bax onto the membrane, leaves it free
# Two step process... should this be comparably fast to tBid binding?
Rule('Bax_to_BaxbraneCplx',
     Bax(b = None, m1 = None, m2 = None, state='I') ** soln +
     Bax(b = 1, m1 = None, m2 = None, state='I') ** lipo_brane % 
     tBid(b1 = 1, b2 = None s='A') ** lipo_brane <>
     Bax(b = 2, m1 = None, m2 = None, state='I') ** lipo_brane %
     Bax(b = 1, m1 = None, m2 = None, state='I') ** lipo_brane % 
     tBid(b = 1, s='A') ** lipo_brane,
     kbaxbranef, kbaxbraner)
Rule('BaxbraneCplx_to_Baxbrane',
     Bax(b = 2, m1 = None, m2 = None, state='I') ** lipo_brane %
     Bax(b = 1  m1 = None, m2 = None, state='I') ** lipo_brane % 
     tBid(b1 = 1, b2 = 2, s='A') ** lipo_brane >>
     Bax(b = None, m1 = None, m2 = None, state='A') ** lipo_brane +
     Bax(b = 1 m1 = None, m2 = None, state='I') ** lipo_brane % 
     tBid(b1 = 1, b2 = 2, s='A') ** lipo_brane,
     kbaxbranec)
     
# Active Bax dimerize
Rule('Baxbrane_dimer',
     Bax(b = None, m1 = None, m2 = None, state='A') ** lipo_brane +
     Bax(b = None, m1 = None, m2 = None, state='A') ** lipo_brane <>
     Bax(b = None, m1 = 1, m2 = None, state='A') ** lipo_brane +
     Bax(b = None, m1 = 1, m2 = None, state='A') ** lipo_brane,
     kdimf, kdimr)
     
# Active Bax trimerize

# Active Bax tetramerize

# Active Bax pentamerize

# Active Bax hexamerize

# Active Bax heptamerize

# CytC mimic released to the cytosol
