from pysb import *

# Basic model aiming to build a model based on the data by 
# Andrews, Chen, and Kinnally

Model()

# Define compartments for this model
Compartment('soln', dimension = 3, size = XX, parent = None)
Compartment('lipo_brane', dimension = 2, size = XX, parent = soln)
Compartment('liposome', dimension = 3, size = XX, parent = lipo_brane)


Monomer('tBid', ['b', 's'], {'s':['act', 'inact']})
Monomer('Bax', ['b', 'm1', 'm2', 's'], {'s':['act', 'inact']})
Monomer('Cc-m', ['b'])

# tBid diffuses to the liposome membrane

# Bax binds to a tBid on the membrane, Bax becomes active

# Active Bax brings another Bax onto the membrane

# Active Bax dimerize

# Active Bax trimerize

# Active Bax tetramerize

# Active Bax pentamerize

# Active Bax hexamerize

# Active Bax heptamerize

# CytC mimic released to the cytosol
