from pysb import *
# Parameter section

# EARM 1.0 parameters
# Reaction rates
parameter_dict = {
    #--------------------
    # EARM 2.0 bcl2 module parameters
    #--------------------
    # Bax transport to mitochondria
    'BAX_trans':  [ Parameter('kbaxCbaxMf', 0.01), 
                    Parameter('kbaxCbaxMr', 0.01)],
    # Bclxl translocation
    'BCLXL_trans':[ Parameter('kbclxlCbclxlMf', 0.01),
                    Parameter('kbclxlCbclxlMr', 0.01)],
    #Bax pore assembly
    'BAX_PORE':   [[Parameter('kbaxdimf',  2.040816e-04),     #1.0e-6/v**2
                    Parameter('kbaxdimr',  1e-3        )],
                   [Parameter('kbaxtrimf', 2.040816e-04),
                    Parameter('kbaxtrimr', 1e-3        )],
                   [Parameter('kbaxtetf',  2.040816e-04),
                    Parameter('kbaxtetr',  1e-3        )]],
    #Bak pore assembly
    'BAK_PORE':   [[Parameter('kbakdimf',  2.040816e-04),
                    Parameter('kbakdimr',  1e-3        )],
                   [Parameter('kbaktrimf', 2.040816e-04),
                    Parameter('kbaktrimr', 1e-3        )],
                   [Parameter('kbaktetf',  2.040816e-04),
                    Parameter('kbaktetr',  1e-3        )]],
    # Inhibitions of Bax/Bak by Bcl2/BclxL/Mcl1
    # These are used in the simple_bind_table function which expects
    # row-major order (if you don't know what this means google it)
    'BID_BAX_BAK_inh':[[Parameter('bidbclxlf', 1.428571e-05),   # 1.0e-6/v
                        Parameter('bidbclxlr', 1.0e-3      )],
                       [Parameter('baxbclxlf', 1.428571e-05),
                        Parameter('baxbclxlr', 1.0e-3      )],
                       [Parameter('bakbclxlf', 1.428571e-05),
                        Parameter('bakbclxlr', 1.0e-3      )],
                       [Parameter('bidmcl1f',  1.428571e-05),
                        Parameter('bidmcl1r',  1.0e-3      )],
                       [Parameter('bakmcl1f',  1.428571e-05),
                        Parameter('bakmcl1r',  1.0e-3      )]],
    # Sensitizers of Bcl2/BclxL/Mcl1 by Bad/NOXA
    'BCLs_sens':      [[Parameter('kbadbclxlf',  1.428571e-05),
                        Parameter('kbadbclxlr',  1.0e-3      )],
                       [Parameter('knoxamcl1f',  1.428571e-05),
                        Parameter('knoxamcl1r',  1.0e-3      )]],
     # CytoC transport by Bax
    'BAX_CYTC':   [[ Parameter('kbaxcytocMCf',  2.857143e-05),   #same as Bax-SMAC
                     Parameter('kbaxcytocMCr',  1.0e-03     ),   #same as Bax-SMAC
                     Parameter('kbaxcytocMCc',  1.0e01      )]], #same as Bax-SMAC
    # Smac transport by Bax
    'BAX_SMAC':   [[ Parameter('kbaxsmacCAf',   2.857143e-05),  #change to MC for consistency (mito to cyto)
                     Parameter('kbaxsmacCAr',   1.0e-03     ),
                     Parameter('kbaxsmacCAc',   1.0e01      )]],
    # CytoC transport activation by Bak **
    'BAK_CYTC':   [[ Parameter('kbakcytocMCf',  2.857143e-05),
                     Parameter('kbakcytocMCr',  1.0e-03     ),
                     Parameter('kbakcytocMCc',  1.0e01      )]],
    # Smac transport by Bak  **
    'BAK_SMAC':   [[ Parameter('kbaksmacCAf',  2.857143e-05),   #same as Bak-SMAC
                     Parameter('kbaksmacCAr',  1.0e-03     ),   #same as Bak-SMAC
                     Parameter('kbaksmacCAc',  1.0e01      )]], #same as Bak-SMAC
    # CytoC activation
    'CYTOC_ACT':  [ Parameter('kcytocCcytoAf', 0.01),  #same as SMAC_ACT
                    Parameter('kcytocCcytoAr', 0.01)], #same as SMAC_ACT
    # Smac activation
    'SMAC_ACT':  [ Parameter('ksmacCsmacAf',   0.01),
                   Parameter('ksmacCsmacAr',   0.01)],
   
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
    'C8_BID':     [ Parameter('kc8bidf', 1.0e-7),
                    Parameter('kc8bidr', 1.0e-3),
                    Parameter('kc8bidc', 1.0)],
    # DISC inhibition by FLIP
    'DISC_FLIP':  [ Parameter('kflipdiscf', 1e-06), 
                    Parameter('kflipdiscr', 1e-03)],
    # C8 inhibition by BAR
    'BAR_C8':     [ Parameter('kbarc8f', 1e-06), 
                    Parameter('kbarc8r', 1e-03)],

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
    'INIT_AMTS':    [ Parameter('L_0'     ,  3000), # 3000 Ligand corresponding to 50 ng/ml SuperKiller TRAIL
                      Parameter('R_0'     ,   200),   # TRAIL receptor 
                      Parameter('flip_0'  , 1.0e2),  # Flip
                      Parameter('C8_0'    , 2.0e4),   # procaspase-8 
                      Parameter('BAR_0'   , 1.0e3),  # Bifunctional apoptosis regulator
                      Parameter('Bid_0'   , 4.0e4),  # Bid
                      Parameter('Bax_0'   , 0.0), # 0.8e5 Bax
                      Parameter('Bax_BclxL_0', 0.8e5), # bax + bclxl
                      Parameter('Bak_0'   , 0.0), #  0.2e5 Bak
                      Parameter('Bak_Mcl1_0', 0.2e5), # bak + mcl1
                      Parameter('BclxL_0' , 2.0e4),  # 2.0e4 cytosolic BclxL
                      Parameter('Mcl1_0'  , 2.0e4),   # 2.0e4 mitochondrial Mcl1  
                      Parameter('Bad_0'   , 1.0e3),  # 1.0e3 Bad
                      Parameter('NOXA_0'  , 1.0e3),  # 1.0e3 NOXA
                      Parameter('CytoC_0' , 5.0e5),  # cytochrome c
                      Parameter('Smac_0'  , 1.0e5),  # Smac    
                      Parameter('Apaf_0'  , 1.0e5),  # Apaf-1
                      Parameter('C3_0'    , 1.0e4),  # 1.0e4 procaspase-3 (pro-C3)
                      Parameter('C6_0'    , 1.0e4),  # procaspase-6 (pro-C6)  
                      Parameter('C9_0'    , 1.0e5),  # procaspase-9 (pro-C9)
                      Parameter('XIAP_0'  , 1.0e5),  # X-linked inhibitor of apoptosis protein  
                      Parameter('PARP_0'  , 1.0e6),  # C3* substrate
                      ]
    }



