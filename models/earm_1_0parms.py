from pysb import *
# Parameter section

# Special parameters
transloc = .01; # rate of transloc bw cytosol and mitochondria
v = .07; # mitochondria compartment volume/cell volume

# Reaction rates
Parameter('klrf', 4e-07) #ligand - receptor forward, reverse, catalytic
Parameter('klrr', 1e-03)
Parameter('klrc', 1e-05)
Parameter('kflipdiscf', 1e-06) #flip-DISC binding
Parameter('kflipdiscr', 1e-03)
Parameter('kdiscc8f', 1e-06) # C8 activation via DISC
Parameter('kdiscc8r', 1e-03)
Parameter('kdiscc8c', 1e+00)
Parameter('kbarc8f', 1e-06) # BAR C8
Parameter('kbarc8r', 1e-03)
Parameter('kc8c3f', 1e-07) # C3 activation
Parameter('kc8c3r', 1e-03)
Parameter('kc8c3c', 1e+00)
Parameter('kc3c6f', 1e-06) # C6 activation
Parameter('kc3c6r', 1e-03)
Parameter('kc3c6c', 1e+00)
Parameter('kc6c8f', 3e-08) # C8 activation via C6
Parameter('kc6c8r', 1e-03)
Parameter('kc6c8c', 1e+00)
Parameter('kxiapc3f', 2e-06) # C3 ubiquitination
Parameter('kxiapc3r', 1e-03)
Parameter('kxiapc3c', 1e-01)
Parameter('kc3parpf', 1e-06) # PARP cleavage by C3
Parameter('kc3parpr', 1e-02)
Parameter('kc3parpc', 1e+00)
Parameter('kc8bidf', 1e-07) # Bid cleavage by C8
Parameter('kc8bidr', 1e-03)
Parameter('kc8bidc', 1e+00)
Parameter('kbidbcl2f', 1e-06) # Bid inhibited by Bcl2
Parameter('kbidbcl2r', 1e-03)
Parameter('kbidbaxf', 1e-07) # Bax activation via Bid
Parameter('kbidbaxr', 1e-03)
Parameter('kbidbaxc', 1e+00)
Parameter('kbaxcbaxmf', transloc) # Bax translocation
Parameter('kbaxcbaxmr', transloc)
Parameter('kbaxMbcl2Mf', 1e-06/v) # Bax inhibition in mito
Parameter('kbaxMbcl2Mr', 1e-03)
Parameter('kbaxdimf', 1e-06/v*2) # Bax dimerization
Parameter('kbaxdimr', 1e-03)

# Initial concentrations
# Non-zero initial conditions (in molecules per cell):
Parameter('L_0'        , 3000); # baseline level of ligand for most experiments (corresponding to 50 ng/ml SuperKiller TRAIL)
Parameter('pR_0'       , 200);  # TRAIL receptor (for experiments not involving siRNA)
Parameter('flip_0'     , 1e2);  # Flip
Parameter('pC8_0'      , 2e4);  # procaspase-8 (pro-C8)
Parameter('BAR_0'      , 1e3);  # Bifunctional apoptosis regulator
Parameter('pC3_0'      , 1e4);  # procaspase-3 (pro-C3)
Parameter('pC6_0'      , 1e4);  # procaspase-6 (pro-C6)  
Parameter('XIAP_0'     , 1e5);  # X-linked inhibitor of apoptosis protein  
Parameter('PARP_0'     , 1e6);  # C3* substrate
Parameter('Bid_0'      , 4e4);  # Bid
Parameter('Bcl2c_0'    , 2e4);  # cytosolic Bcl-2
Parameter('Bax_0'      , 1e5);  # Bax
Parameter('Bcl2_0'     , 2e4);  # mitochondrial Bcl-2  
Parameter('Mito_0'     , 5e5);  # mitochondrial binding sites for activated Bax
Parameter('mCytoC_0'   , 5e5);  # cytochrome c
Parameter('mSmac_0'    , 1e5);  # Smac    
Parameter('pC9_0'      , 1e5);  # procaspase-9 (pro-C9)
Parameter('Apaf_0'     , 1e5);  # Apaf-1

