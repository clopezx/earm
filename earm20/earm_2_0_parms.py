from pysb import *
# Parameter section
#
#

# Old parameters parameters
#v = .07; # mitochondria compartment volume/cell volume
#tranloc = .01

# EARM 1.0 parameters
# Reaction rates
parameter_dict = {
    #--------------------
    # EARM 2.0 bcl2 module parameters, allowed to vary O(3) during annealing
    #--------------------
    # Bid transport to mitochondria
    'BID_trans':  [ Parameter('kbidTbidMf',  .01), #Lovell, Andrews Cell 2008 says 1e12?
                    Parameter('kbidTbidMr',  .01)],
    # Bad transport to mitochondria
    'BAD_trans':  [ Parameter('kbadCbadMf', .01),
                    Parameter('kbadCbadMr', .01)],
    # Bax transport to mitochondria
    'BAX_trans':  [ Parameter('kbaxCbaxMf',  .01), 
                    Parameter('kbaxCbaxMr',  .01)],
    # Bcl2 translocation
    'BCL2_trans': [ Parameter('kbcl2Cbcl2Mf', .01),
                    Parameter('kbcl2Cbcl2Mr', .01)],
    # Bclxl translocation
    'BCLXL_trans':[ Parameter('kbclxlCbclxlMf', 0.01),
                    Parameter('kbclxlCbclxlMr', 0.01)],
    # Bax activation by Bid
    'BID_BAX':    [ Parameter('kbidbaxf', 1.0e-07),
                    Parameter('kbidbaxr', 1.0e-03),
                    Parameter('kbidbaxc', 1.0    )],
    # Bak activation via Bid ***CHECK VALUES
    'BID_BAK':    [ Parameter('kbidbakf', 1.0e-07),
                    Parameter('kbidbakr', 1.0e-03),
                    Parameter('kbidbakc', 1.0    )],
    # Bax activation by Bid
    'BAX_BAX':    [ Parameter('kbaxbaxf', 1.0e-07),
                    Parameter('kbaxbaxr', 1.0e-03),
                    Parameter('kbaxbaxc', 1.0    )],
    # Bak activation via Bid ***CHECK VALUES
    'BAK_BAK':    [ Parameter('kbakbakf', 1.0e-07),
                    Parameter('kbakbakr', 1.0e-03),
                    Parameter('kbakbakc', 1.0    )],
    #Bax pore assembly
    'BAX_PORE':   [[Parameter('kbaxdimf',  2.040816e-04), #1.0e-6/v**2
                    Parameter('kbaxdimr',  1.0e-3      )],
                   [Parameter('kbaxtrimf', 2.040816e-04),
                    Parameter('kbaxtrimr', 1.0e-3      )],
                   [Parameter('kbaxtetf',  2.040816e-04),
                    Parameter('kbaxtetr',  1.0e-3      )]],
    #Bak pore assembly
    'BAK_PORE':   [[Parameter('kbakdimf',  2.040816e-04),
                    Parameter('kbakdimr',  1.0e-3      )],
                   [Parameter('kbaktrimf', 2.040816e-04),
                    Parameter('kbaktrimr', 1.0e-3      )],
                   [Parameter('kbaktetf',  2.040816e-04),
                    Parameter('kbaktetr',  1.0e-3      )]],
    #BclxL recruitment by and inhibition of Bid
    'Bid_BclxL_RA': [Parameter('kbidbclxl_RAf', 2.040816e-04),
                     Parameter('kbidbclxl_RAr', 1.0e-3      ),
                     Parameter('kbidbclxl_RAc', 1.0         )],
    #BclxL recruitment by and inhibition of Bax
    'Bax_BclxL_RA': [Parameter('kbaxbclxl_RAf', 2.040816e-04),
                     Parameter('kbaxbclxl_RAr', 1.0e-3      ),
                     Parameter('kbaxbclxl_RAc', 1.0         )],
    # Inhibitions of Bax/Bak by Bcl2/BclxL/Mcl1
    # These are used in the simple_bind_table function which expects
    # row-major order (if you don't know what this means google it)
    'BID_BAX_BAK_inh':[[Parameter('kbidbcl2f', 1.428571e-05), #1.0e-06/v
                        Parameter('kbidbcl2r', 1.0e-3      )],
                       [Parameter('baxbcl2f',  1.428571e-05), #
                        Parameter('baxbcl2r',  1.0e-3      )],
                       [Parameter('bidbclxlf', 1.428571e-05),#
                        Parameter('bidbclxlr', 1.0e-3      )],
                       [Parameter('baxbclxlf', 1.428571e-05),#
                        Parameter('baxbclxlr', 1.0e-3      )],
                       [Parameter('bakbclxlf', 1.428571e-05),#
                        Parameter('bakbclxlr', 1.0e-3      )],
                       [Parameter('bidmcl1f',  1.428571e-05),#
                        Parameter('bidmcl1r',  1.0e-3      )],
                       [Parameter('bakmcl1f',  1.428571e-05),#
                        Parameter('bakmcl1r',  1.0e-3      )]],

    # Sensitizers of Bcl2/BclxL/Mcl1 by Bad/NOXA
    'BCLs_sens':      [[Parameter('kbadbcl2f',  1.428571e-05),
                        Parameter('kbadbcl2r',  1.0e-3      )],
                       [Parameter('kbadbclxlf', 1.428571e-05),
                        Parameter('kbadbclxlr', 1.0e-3      )],
                       [Parameter('knoxamcl1f', 1.428571e-05),
                        Parameter('knoxamcl1r', 1.0e-3      )]],
    
     # CytoC transport by Bax
    'BAX_CYTC':   [[ Parameter('kbaxcytocMCf', 2.857143e-05),   #same as Bax-SMAC
                     Parameter('kbaxcytocMCr', 1.0e-03     ),   #same as Bax-SMAC
                     Parameter('kbaxcytocMCc', 1.0         )]], #same as Bax-SMAC
    # Smac transport by Bax
    'BAX_SMAC':   [[ Parameter('kbaxsmacCAf',  2.857143e-05),
                     Parameter('kbaxsmacCAr',  1.0e-03     ),
                     Parameter('kbaxsmacCAc',  1.0         )]],
    # CytoC transport activation by Bak 
    'BAK_CYTC':   [[ Parameter('kbakcytocMCf', 2.857143e-05),
                     Parameter('kbakcytocMCr', 1.0e-03     ),
                     Parameter('kbakcytocMCc', 1.0         )]],
    # Smac transport by Bak  **
    'BAK_SMAC':   [[ Parameter('kbaksmacCAf',  2.857143e-05),   #same as Bak-SMAC
                     Parameter('kbaksmacCAr',  1.0e-03     ),   #same as Bak-SMAC
                     Parameter('kbaksmacCAc',  1.0         )]], #same as Bak-SMAC
    # CytoC activation
    'CYTOC_ACT':  [ Parameter('kcytocCcytoAf', 0.01), #transloc
                    Parameter('kcytocCcytoAr', 0.01)], 
    # Smac activation
    'SMAC_ACT':  [ Parameter('ksmacCsmacAf', .01), #transloc
                   Parameter('ksmacCsmacAr', .01)],

    #---------------------------
    # EARM 1.0 legacy parameters
    #---------------------------

    # rec_to_bid module parameters
    'L_R_DISC':   [ Parameter('klrf', 4.0e-07),       
                    Parameter('klrr', 1.0e-03),
                    Parameter('klrc', 1.0e-05)],
    # C8 activation via DISC,
    'DISC_C8':    [ Parameter('kdiscc8f', 1.0e-06),   
                    Parameter('kdiscc8r', 1.0e-03),
                    Parameter('kdiscc8c', 1.0    )],
    # Bid cleavage by C8
    'C8_BID':     [ Parameter('kc8bidf', 1.0e-07),
                    Parameter('kc8bidr', 1.0e-03),
                    Parameter('kc8bidc', 1.0    )],
    # DISC inhibition by FLIP
    'DISC_FLIP':  [ Parameter('kflipdiscf', 1.0e-06), 
                    Parameter('kflipdiscr', 1.0e-03)],
    # C8 inhibition by BAR
    'BAR_C8':     [ Parameter('kbarc8f', 1.0e-06), 
                    Parameter('kbarc8r', 1.0e-03)],

    #---------------------------
    # pore_to_parp module parameters
    #---------------------------
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
    'SMAC_XIAP':  [ Parameter('ksmacxiapf', 7.0e-06),
                    Parameter('ksmacxiapr', 1.0e-03)],
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
                      Parameter('R_0'     ,   200), # TRAIL receptor 
                      Parameter('flip_0'  , 1.0e2), # Flip
                      Parameter('C8_0'    , 2.0e4), # procaspase-8 
                      Parameter('BAR_0'   , 1.0e3), # Bifunctional apoptosis regulator
                      Parameter('Bid_0'   , 4.0e4), # Bid
                      Parameter('Bax_0'   , 0.8e5), # Bax
                      Parameter('Bak_0'   , 0.2e5), # Bak
                      Parameter('Bcl2_0'  , 2.0e4), # cytosolic Bcl2
                      Parameter('BclxL_0' , 2.0e4), # cytosolic BclxL
                      Parameter('Mcl1_0'  , 2.0e4), # mitochondrial Mcl1  
                      Parameter('Bad_0'   , 1.0e3), # Bad
                      Parameter('NOXA_0'  , 1.0e3), # NOXA
                      Parameter('CytoC_0' , 5.0e5), # cytochrome c
                      Parameter('Smac_0'  , 1.0e5), # Smac    
                      Parameter('Apaf_0'  , 1.0e5), # Apaf-1
                      Parameter('C3_0'    , 1.0e4), # procaspase-3 (pro-C3)
                      Parameter('C6_0'    , 1.0e4), # procaspase-6 (pro-C6)  
                      Parameter('C9_0'    , 1.0e5), # procaspase-9 (pro-C9)
                      Parameter('XIAP_0'  , 1.0e5), # X-linked inhibitor of apoptosis protein  
                      Parameter('PARP_0'  , 1.0e6), # C3* substrate
                      ]
    }


