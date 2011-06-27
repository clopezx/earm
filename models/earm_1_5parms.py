from pysb import *
# Parameter section

# Special parameters
transloc = .01; # rate of transloc bw cytosol and mitochondria
v = .07; # mitochondria compartment volume/cell volume

# EARM 1.0 parameters
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
Parameter('kbaxCbaxMf', transloc) # Bax translocation
Parameter('kbaxCbaxMr', transloc)
Parameter('kbidCbidMf', transloc) # Bid translocation
Parameter('kbidCbidMr', transloc)
Parameter('kbaxMbcl2Mf', 1e-06/v) # Bax inhibition in mito
Parameter('kbaxMbcl2Mr', 1e-03)
Parameter('kbaxdimf', 1e-06/v*2) # Bax dimerization
Parameter('kbaxdimr', 1e-03)
Parameter('kbax2Mbcl2Mf', 1e-06/v) # Bax2 inibition in mito
Parameter('kbax2Mbcl2Mr', 1e-03)
Parameter('kbaxtetf', 1e-06/v*2) # Bax2 dimerization (tetramer formation)
Parameter('kbaxtetr', 1e-03)
Parameter('kbax4Mbcl2Mf', 1e-06/v) # Bax4 inhibition
Parameter('kbax4Mbcl2Mr', 1e-03)
Parameter('kbax4poref', 1e-06/v) # Bax4 + MitoP to Pore
Parameter('kbax4porer', 1e-03)
Parameter('kbax4porec', 1e+00)
Parameter('kmitopcytocMf', 2e-06/v) # CytoC activation
Parameter('kmitopcytocMr', 1e-03)
Parameter('kmitopcytocMc', 1e+01)
Parameter('kmitopsmacMf', 2e-06/v) # Smac activation
Parameter('kmitopsmacMr', 1e-03)
Parameter('kmitopsmacMc', 1e+01)
Parameter('kcytocMcytocCf', transloc) # CytoC translocation
Parameter('kcytocMcytocCr', transloc)
Parameter('kcytocCapaff', 5e-07) # Apaf activation
Parameter('kcytocCapafr', 1e-03)
Parameter('kcytocCapafc', 1e+00)
Parameter('kapafc9f', 5e-08) # Apop formation
Parameter('kapafc9r', 1e-03)
Parameter('kapopc3f', 5e-09) # C3 activation via Apop
Parameter('kapopc3r', 1e-03)
Parameter('kapopc3c', 1e+00)
Parameter('ksmacMsmacCf', transloc) # Smac translocation
Parameter('ksmacMsmacCr', transloc)
Parameter('kapopxiapf', 2e-06) # Apop inhibition by XIAP
Parameter('kapopxiapr', 1e-03)
Parameter('ksmacxiapf', 7e-06) # XIAP inhibition by Smac
Parameter('ksmacxiapr', 1e-03)

# EARM 1.5 parameters
Parameter('kbcl2Cbcl2Mf', transloc) # Bcl2 translocation
Parameter('kbcl2Cbcl2Mr', transloc)
Parameter('kbclxlCbclxlMf', transloc) # Bclxl translocation
Parameter('kbclxlCbclxlMr', transloc)
Parameter('baxbcl2f', 3.33) # Inhibitions of Bax/Bak by Bcl2/BclxL/Mcl1
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
Parameter('badbcl2f', 3.33) # Sensitization of Bcl2/BclxL/Mcl1 by Bad/NOXA
Parameter('badbcl2r', 3.33)
Parameter('badbclxlf', 3.33) 
Parameter('badbclxlr', 3.33)
Parameter('badmcl1f', 3.33)
Parameter('badmcl1r', 3.33)
Parameter('noxabcl2f', 3.33)
Parameter('noxabcl2r', 3.33)
Parameter('noxabclxlf', 3.33)
Parameter('noxabclxlr', 3.33)
Parameter('noxamcl1f', 3.33)
Parameter('noxamcl1r', 3.33)


# Initial amounts
# Non-zero initial conditions (in molecules per cell):
Parameter('L_0'        , 3000); # baseline level of ligand for most experiments (corresponding to 50 ng/ml SuperKiller TRAIL)
Parameter('R_0'       , 200);  # TRAIL receptor (for experiments not involving siRNA)
Parameter('flip_0'     , 1e2);  # Flip
Parameter('C8_0'      , 2e4);  # procaspase-8 (pro-C8)
Parameter('BAR_0'      , 1e3);  # Bifunctional apoptosis regulator
Parameter('Bid_0'      , 4e4);  # Bid
Parameter('Bax_0'      , 1e5);  # Bax
Parameter('Bak_0'      , 1e0);  # Bax
Parameter('Bcl2_0'    , 2e4);  # cytosolic Bcl2
Parameter('BclxL_0'    , 2e4);  # cytosolic BclxL
Parameter('Mcl1_0', 2e4);  # mitochondrial Mcl1  
Parameter('Bad_0'      , 1e3);  # Bad
Parameter('NOXA_0'      , 1e3);  # NOXA
Parameter('CytoC_0'   , 5e5);  # cytochrome c
Parameter('Smac_0'    , 1e5);  # Smac    
Parameter('Apaf_0'     , 1e5);  # Apaf-1
Parameter('C3_0'      , 1e4);  # procaspase-3 (pro-C3)
Parameter('C6_0'      , 1e4);  # procaspase-6 (pro-C6)  
Parameter('C9_0'      , 1e5);  # procaspase-9 (pro-C9)
Parameter('XIAP_0'     , 1e5);  # X-linked inhibitor of apoptosis protein  
Parameter('PARP_0'     , 1e6);  # C3* substrate
