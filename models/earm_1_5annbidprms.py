from pysb import *
# Parameter section
#
# 8/3/2011 
# The parameters for ['kbidCbidMf', 'kbidCbidMr', 'kbidbcl2f', 'kbidbcl2r', 'kc8bidf', 'kc8bidr', 'kc8bidc']
# were optimized using simulated annealing to match the output of tBid in earm_1_5 to the 
# original output of earm_1_0. I used a sigma of 5% (i.e. .05) for the objective function.
#

# Special parameters
v = .07; # mitochondria compartment volume/cell volume

# EARM 1.0 parameters
# Reaction rates
parameter_dict = {
    #--------------------
    # EARM 1.5 bcl2 module parameters
    #--------------------
    # Bid transport to mitochondria
    'BID_trans':  [ Parameter('kbidCbidMf', 4.91481326e+00),
                    Parameter('kbidCbidMr', 1.40411075e-03)],
    # Bax transport to mitochondria
    'BAX_trans':  [ Parameter('kbaxCbaxMf', .01), 
                    Parameter('kbaxCbaxMr', .01)],
    # Bcl2 translocation
    'BCL2_trans': [ Parameter('kbcl2Cbcl2Mf', .01),
                    Parameter('kbcl2Cbcl2Mr', .01)],
    # Bclxl translocation
    'BCLXL_trans':[ Parameter('kbclxlCbclxlMf', .01),
                    Parameter('kbclxlCbclxlMr', .01)],
    # Bax activation by Bid
    'BID_BAX':    [ Parameter('kbidbaxf', 1e-07),
                    Parameter('kbidbaxr', 1e-03),
                    Parameter('kbidbaxc', 1e+00)],
    # Bak activation via Bid ***CHECK VALUES
    'BID_BAK':    [ Parameter('kbidbakf', 1e-07),
                    Parameter('kbidbakr', 1e-03),
                    Parameter('kbidbakc', 1e+00)],
    #Bax pore assembly
    'BAX_PORE':   [[Parameter('kbaxdimf', 1e-06/v*2),
                    Parameter('kbaxdimr', 1e-03)],
                   [Parameter('kbaxtrimf', 1e-06/v*2),
                    Parameter('kbaxtrimr', 1e-03)],
                   [Parameter('kbaxtetf', 1e-06/v*2),
                    Parameter('kbaxtetr', 1e-03)]],
    'BAK_PORE':   [[Parameter('kbakdimf', 1e-06/v*2),
                    Parameter('kbakdimr', 1e-03)],
                   [Parameter('kbaktrimf', 1e-06/v*2),
                    Parameter('kbaktrimr', 1e-03)],
                   [Parameter('kbaktetf', 1e-06/v*2),
                    Parameter('kbaktetr', 1e-03)]],
    # Inhibitions of Bax/Bak by Bcl2/BclxL/Mcl1
    # These are used in the simple_bind_table function which expects
    # row-major order (if you don't know what this means google it)
    'BID_BAX_BAK_inh':[[Parameter('kbidbcl2f', 1.45280010e-06),
                        Parameter('kbidbcl2r', 1.70789566e-03)],
                       [Parameter('baxbcl2f', 1e-05/v), #
                        Parameter('baxbcl2r', 1e-03  )],
                       [Parameter('baxbclxlf', 1e-05/v),#
                        Parameter('baxbclxlr', 1e-03  )],
                       [Parameter('bakbclxlf',1e-05/v),#
                        Parameter('bakbclxlr',1e-03  )],
                       [Parameter('bakmcl1f', 1e-05/v),#
                        Parameter('bakmcl1r', 1e-03  )]],
    # Sensitizers of Bcl2/BclxL/Mcl1 by Bad/NOXA
    'BCLs_sens':      [[Parameter('kbadbcl2f',  1e-09/v),
                        Parameter('kbadbcl2r',  1e-03  )],
                       [Parameter('kbadbclxlf', 1e-09/v),
                        Parameter('kbadbclxlr', 1e-03  )],
                       [Parameter('knoxabcl2f', 1e-09/v),
                        Parameter('knoxabcl2r', 1e-03  )],
                       [Parameter('knoxamcl1f', 1e-09/v),
                        Parameter('knoxamcl1r', 1e-03  )]],
    
    #---------------------------
    # EARM 1.0 legacy parameters
    #---------------------------
    # rec_to_bid module parameters
    'L_R_DISC':   [ Parameter('klrf', 4e-07),       
                    Parameter('klrr', 1e-03),
                    Parameter('klrc', 1e-05)],
    # C8 activation via DISC,
    'DISC_C8':    [ Parameter('kdiscc8f', 1e-06),   
                    Parameter('kdiscc8r', 1e-03),
                    Parameter('kdiscc8c', 1e+00)],
    # Bid cleavage by C8
    'C8_BID':     [ Parameter('kc8bidf', 9.98683170e-08),
                    Parameter('kc8bidr', 1.72707515e-03),
                    Parameter('kc8bidc', 1.40842965e+00)],
    # DISC inhibition by FLIP
    'DISC_FLIP':  [ Parameter('kflipdiscf', 1e-06), 
                    Parameter('kflipdiscr', 1e-03)],
    # C8 inhibition by BAR
    'BAR_C8':     [ Parameter('kbarc8f', 1e-06), 
                    Parameter('kbarc8r', 1e-03)],

    #---------------------------
    # pore_to_parp module parameters
    #---------------------------
    # CytoC transport/activation by Bax
    'BAX_CYTC':   [[ Parameter('kbaxcytocMCf', 2e-06/v),
                     Parameter('kbaxcytocMCr', 1e-03),
                     Parameter('kbaxcytocMCc', 1e+01)]],
    # Smac transport/activation by Bax
    'BAX_SMAC':   [[ Parameter('kbaxsmacCAf', 2e-06/v),
                     Parameter('kbaxsmacCAr', 1e-03),
                     Parameter('kbaxsmacCAc', 1e+01)]],
    # CytoC transport/activation by Bak **
    'BAK_CYTC':   [[ Parameter('kbakcytocMCf', 2e-06/v),
                     Parameter('kbakcytocMCr', 1e-03),
                     Parameter('kbakcytocMCc', 1e+01)]],
    # Smac transport/activation by Bak  **
    'BAK_SMAC':   [[ Parameter('kbaksmacCAf', 2e-06/v),
                     Parameter('kbaksmacCAr', 1e-03),
                     Parameter('kbaksmacCAc', 1e+01)]],
    # Apaf activation by CytC
    'APAF_CYTC':  [ Parameter('kcytocCapaff', 5e-07),
                    Parameter('kcytocCapafr', 1e-03),
                    Parameter('kcytocCapafc', 1e+00)],
    # Apop formation by Apaf + C9
    'APOP_C9:APAF':[ Parameter('kapafc9f', 5e-08),
                     Parameter('kapafc9r', 1e-03)],
     # C3 activation by Apop
    'APOP_C3':    [ Parameter('kapopc3f', 5e-09),
                    Parameter('kapopc3r', 1e-03),
                    Parameter('kapopc3c', 1e+00)],
    # Apop inhibition by XIAP
    'APOP_XIAP':  [ Parameter('kapopxiapf', 2e-06),
                    Parameter('kapopxiapr', 1e-03)],
    # XIAP inhibition by Smac
    'SMAC_XIAP':  [ Parameter('ksmacxiapf', 7e-06),
                    Parameter('ksmacxiapr', 1e-03)],
    # C3 activation by C8
    'C3_C8':      [ Parameter('kc8c3f', 1e-07), 
                    Parameter('kc8c3r', 1e-03),
                    Parameter('kc8c3c', 1e+00)],
    # C6 activation by C3
    'C6_C3':      [ Parameter('kc3c6f', 1e-06), 
                    Parameter('kc3c6r', 1e-03),
                    Parameter('kc3c6c', 1e+00)],
    # C3 ubiquination tag by XIAP
    'C3_XIAP':    [ Parameter('kxiapc3f', 2e-06),
                    Parameter('kxiapc3r', 1e-03),
                    Parameter('kxiapc3c', 1e-01)],
    # PARP cleavage by C3
    'PARP_C3':    [ Parameter('kc3parpf', 1e-06),
                    Parameter('kc3parpr', 1e-02),
                    Parameter('kc3parpc', 1e+00)],
    # C8 activation by C6
    'C8_C6':      [ Parameter('kc6c8f', 3e-08), 
                    Parameter('kc6c8r', 1e-03),
                    Parameter('kc6c8c', 1e+00)],

    #---------------------------------
    # EARM 1.0 HeLa initial conditions
    # Non-zero initial conditions (in molecules per cell):
    #---------------------------------
    'INIT_AMTS':    [ Parameter('L_0'     ,  3000), # Ligand corresponding to 50 ng/ml SuperKiller TRAIL
                      Parameter('R_0'     ,   200),   # TRAIL receptor 
                      Parameter('flip_0'  , 1.0e2),  # Flip
                      Parameter('C8_0'    , 2.0e4),   # procaspase-8 
                      Parameter('BAR_0'   , 1.0e3),  # Bifunctional apoptosis regulator
                      Parameter('Bid_0'   , 4.0e4),  # Bid
                      Parameter('Bax_0'   , 0.8e5), # Bax
                      Parameter('Bak_0'   , 0.2e5), # Bak
                      Parameter('Bcl2_0'  , 2.0e4),   # cytosolic Bcl2
                      Parameter('BclxL_0' , 2.0e4),  # cytosolic BclxL
                      Parameter('Mcl1_0'  , 2.0e4),   # mitochondrial Mcl1  
                      Parameter('Bad_0'   , 1.0e3),  # Bad
                      Parameter('NOXA_0'  , 1.0e3),  # NOXA
                      Parameter('CytoC_0' , 5.0e5),  # cytochrome c
                      Parameter('Smac_0'  , 1.0e5),  # Smac    
                      Parameter('Apaf_0'  , 1.0e5),  # Apaf-1
                      Parameter('C3_0'    , 1.0e4),  # procaspase-3 (pro-C3)
                      Parameter('C6_0'    , 1.0e4),  # procaspase-6 (pro-C6)  
                      Parameter('C9_0'    , 1.0e5),  # procaspase-9 (pro-C9)
                      Parameter('XIAP_0'  , 1.0e5),  # X-linked inhibitor of apoptosis protein  
                      Parameter('PARP_0'  , 1.0e6),  # C3* substrate
                      ]
    }



