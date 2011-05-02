from pysb import *

# Basic model aiming to build a model based on the data by 
# Andrews, Chen, and Kinnally

Model()

# Compartments section
Compartment('soln', dimension = 3, size = XX, parent = None)
Compartment('lipo_brane', dimension = 2, size = XX, parent = soln)
Compartment('liposome', dimension = 3, size = XX, parent = lipo_brane)

# Species monomer sections
Monomer('tBid', ['b', 'state'], {'state':['act', 'inact']})
Monomer('Bax', ['b', 'm1', 'm2', 'state'], {'state':['act', 'inact']})
Monomer('Cc-m', ['b'])

Parameter('[NAME]', [VALUE])

# tBid diffuses to the liposome membrane
Rule('tBid_to_brane',
     tBid(b = None, s='act') ** soln <> tBid(b = None, s='act') ** lipo_brane,
     ktbidbranef, ktbidbraner)

# Bax binds to a tBid on the membrane, Bax becomes active
Rule('Bax_to_Bidbrane',
     Bax(b = None, m1 = None, m2 = None, state='inact') ** soln + 
     tBid(b = None, s='act') ** libo_brane <>
     Bax(b = 1, m1 = None, m2 = None, state='act') ** lipo_brane % 
     tBid(b = 1, s='act') ** libo_brane,
     kbaxbranef, baxbraner)

# Active Bax brings another Bax onto the membrane, leaves it free
# FIXME: Make this a complex then products reaction rather than a one-step thing?
Rule('Bax_to_Baxbrane',
     Bax(b = None, m1 = None, m2 = None, state='inact') ** soln +
     Bax(b = 1, m1 = None, m2 = None, state='act') ** lipo_brane % 
     tBid(b = 1, s='act') ** libo_brane <>
     Bax(b = None, m1 = None, m2 = None, state='act') ** lipo_brane +
     Bax(b = 1, m1 = None, m2 = None, state='act') ** lipo_brane % 
     tBid(b = 1, s='act') ** libo_brane,
     kbaxbranef, kbaxbraner)
     
     


# Active Bax dimerize

# Active Bax trimerize

# Active Bax tetramerize

# Active Bax pentamerize

# Active Bax hexamerize

# Active Bax heptamerize

# CytC mimic released to the cytosol
