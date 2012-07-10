from pysb import *
from collections import OrderedDict

# EARM 2.0 embedded parameters
# from simulated annealing fits to experimental data using a sum of chi-sq objective functions
#

#ec_size = 1.0e6             # 1.0e6 um^3 = 1 pL								 
#cytoM_size = 483.6 * .0030  # plasma SA (6.22um radius for a 1e3 um^3 cell) * membrane thickness ~3.0nm	 
#cyto_size = 1.0e3	    # 1.0e3 um^3 --> size of HeLa. Range is 760-2730 um^3 (ref)			 
#mito_size = 70.0	    # mitochondria is ~7% of cytoplasm (ref)					 
#mitoM_size = 82.14 * .0042  # mito SA (2.55um radius) x 'brane thicknes ~4.2 nm JPC-B(2009)113-11p3413   
ec_size = 1.0   
cytoM_size = 1.0
cyto_size = 1.0
mito_size = 1.0
mitoM_size = 1.0

# Compartment parameters
Parameter('ec_size', ec_size)       
Parameter('cytoM_size', cytoM_size) 
Parameter('cyto_size', cyto_size)   
Parameter('mito_size', mito_size)   
Parameter('mitoM_size', mitoM_size) 


# Reaction rates
# Note: these parameters were imported from earm_v2 with implicit volume and must be
#       corrected for a run with explicit volumes
parameter_dict = OrderedDict([
    # Trail to Bid parameters
    # -----------------------
    # rec_to_bid module parameters
    # TRAIL trimer binding to DR monomer- 4, values from Reis + Cool 2011
   ('TT_DR4mono',  [ Parameter('ttdr4mf', 1.040000e-06 / ec_size  ),    
                     Parameter('ttdr4mr', 1.10e-04)]),
    # TT binding to DR dimers- 4
   ('TT_DR4dim',   [ Parameter('ttdr4df', 5.000000e-07 / ec_size  ),    
                     Parameter('ttdr4dr', 1.10e-04)]),
    # TT binding to DR trimers- 4
   ('TT_DR4trim',  [ Parameter('ttdr4tf', 1.040000e-07 / cyto_size  ),     
                     Parameter('ttdr4tr', 1.100000e-04 / cyto_size  )]),
    # TRAIL trimer binding to DR monomer- 5
    # TT binding to DR monomer- 5
   ('TT_DR5mono',  [ Parameter('ttdr5mf', 1.970000e-06 / ec_size  ),    
                     Parameter('ttdr5mr', 0.36e-04)]),
    # TT binding to DR dimers- 5
   ('TT_DR5dim',   [ Parameter('ttdr5df', 2.070000e-06 / ec_size  ),    
                     Parameter('ttdr5dr', 1.10e-04)]),
    # TT binding to DR trimers- 5
   ('TT_DR5trim',  [ Parameter('ttdr5tf', 1.970000e-07 / cyto_size  ),     
                     Parameter('ttdr5tr', 3.600000e-05 / cyto_size  )]),
    # DR4 trimerization
   ('DR4_RINGP',  [[Parameter('kdr4dimf', 2.040816e-04 / cytoM_size  ), 
                    Parameter('kdr4dimr',  1.0e-3      )],
                    [Parameter('kdr4trimf', 2.040816e-04 / cytoM_size  ),
                    Parameter('kdr4trimr', 1.0e-3      )]]),
    # DR5 trimerization
   ('DR5_RINGP',  [[Parameter('kdr5dimf', 2.040816e-04 / cytoM_size  ), 
                    Parameter('kdr5dimr',  1.0e-3      )],
                   [Parameter('kdr5trimf', 2.040816e-04 / cytoM_size  ),
                    Parameter('kdr5trimr', 1.0e-3      )]]),
   # Fadd LDR4C complex
   ('Fadd_LDR4C_C1', [Parameter('faddldr4c1f', 4.000000e-07 / cytoM_size  ),       
                      Parameter('faddldr4c1r', 1.0e-03)]),
   ('Fadd_LDR4C_C2', [Parameter('faddldr4c2f', 4.000000e-07 / cytoM_size  ),       
                      Parameter('faddldr4c2r', 1.0e-03)]),
   ('Fadd_LDR4C_C3', [Parameter('faddldr4c3f', 4.000000e-07 / cytoM_size  ),       
                      Parameter('faddldr4c3r', 1.0e-03)]),
   # Fadd LDR5C complex
   ('Fadd_LDR5C_C1', [Parameter('faddldr5c1f', 4.000000e-07 / cytoM_size  ),       
                      Parameter('faddldr5c1r', 1.0e-03)]),
   ('Fadd_LDR5C_C2', [Parameter('faddldr5c2f', 4.000000e-07 / cytoM_size  ),       
                      Parameter('faddldr5c2r', 1.0e-03)]),
   ('Fadd_LDR5C_C3', [Parameter('faddldr5c3f', 4.000000e-07 / cytoM_size  ),       
                      Parameter('faddldr5c3r', 1.0e-03)]),
    # pC8 FADD binding pC8_fadd_b
   ('pC8_fadd_b', [ Parameter('pc8faddf', 4.000000e-07 / cyto_size  ),       
                    Parameter('pc8faddr', 1.0e-03)]),
    # pC8 activation to C8 within the same FADD trimer
   ('pC8_dim_act_s', [ Parameter('pc8dimacts', 1.0e-3)]),
    # pC8 activation to C8 across separate FADD trimers
   ('pC8_dim_act_o', [ Parameter('pc8dimacto', 1.0e-3)]),
    # Bid cleavage by C8
   ('C8_BID',     [ Parameter('kc8bidf', 1.000000e-07 / cyto_size  ),
                    Parameter('kc8bidr', 1.0e-03),
                    Parameter('kc8bidc', 1.0    )]),
    # DISC inhibition by FLIP
   ('DISC_FLIP',  [ Parameter('kflipdiscf', 1.000000e-06 / cytoM_size  ), 
                    Parameter('kflipdiscr', 1.0e-03)]),
    # C8 inhibition by BAR
   ('BAR_C8',     [ Parameter('kbarc8f', 1.0e-06), 
                    Parameter('kbarc8r', 1.0e-03)]),

    # BID TO MOMP parameters
    # ----------------------
    #
                    
    ( 'BID_trans' ,
      [Parameter('kbidTbidMf', 0.0026601962092430266), Parameter('kbidTbidMr', 6.4294317349802732e-05)] ),
    ( 'BAD_trans' ,
      [Parameter('kbadCbadMf', 0.0029419710017770154), Parameter('kbadCbadMr', 0.0027381917924608314)] ),
    ( 'BAX_trans' ,
      [Parameter('kbaxCbaxMf', 0.0027070779785728404), Parameter('kbaxCbaxMr', 0.004353093828158041)] ),
    ( 'BCL2_trans' ,
      [Parameter('kbcl2Cbcl2Mf', 0.0014737610524655136), Parameter('kbcl2Cbcl2Mr', 0.014906933079643386)] ),
    ( 'BCLXL_trans' ,
      [Parameter('kbclxlCbclxlMf', 0.027466171257673082), Parameter('kbclxlCbclxlMr', 0.0084157689975014327)] ),
    ( 'BID_BAX' ,
      [Parameter('kbidbaxf', 3.895653e-07 / cyto_size  ),
       Parameter('kbidbaxr', 0.00036856169047138091),
       Parameter('kbidbaxc', 1.4664186386875186)] ),
    ( 'BID_BAK' ,
      [Parameter('kbidbakf', 6.155795e-09 / mitoM_size  ),
       Parameter('kbidbakr', 0.0029106994159151453),
       Parameter('kbidbakc', 0.84528236986831273)] ),
    ( 'BAX_BAX' ,
      [Parameter('kbaxbaxf', 3.046826e-10 / cyto_size  ),
       Parameter('kbaxbaxr', 1.7249510431983885e-05),
       Parameter('kbaxbaxc', 2.1191713445212592)] ),
    ( 'BAK_BAK' ,
      [Parameter('kbakbakf', 3.663557e-08 / mitoM_size  ),
       Parameter('kbakbakr', 0.0013849096067022136),
       Parameter('kbakbakc', 0.23567721868317815)] ),
    ( 'BAX_PORE' ,
      [[Parameter('kbaxdimf', 1.987356e-04 / mitoM_size  ), Parameter('kbaxdimr', 0.00051552975969350876)],
       [Parameter('kbaxtrimf', 2.674313e-04 / mitoM_size  ), Parameter('kbaxtrimr', 0.00053921595886761769)],
       [Parameter('kbaxtetf', 5.333138e-05 / mitoM_size  ), Parameter('kbaxtetr', 1.6982589396135685e-05)]] ),
    ( 'BAK_PORE' ,
      [[Parameter('kbakdimf', 4.438236e-04 / mitoM_size  ), Parameter('kbakdimr', 0.0017427126140823314)],
       [Parameter('kbaktrimf', 7.306882e-04 / mitoM_size  ), Parameter('kbaktrimr', 0.00041757049236862946)],
       [Parameter('kbaktetf', 1.585759e-05 / mitoM_size  ), Parameter('kbaktetr', 0.0010660677426411207)]] ),
    ( 'Bid_BclxL_RA' ,
      [Parameter('kbidbclxl_RAf', 6.934216e-04 / cyto_size  ),
       Parameter('kbidbclxl_RAr', 0.00083046693481875049),
       Parameter('kbidbclxl_RAc', 0.69567061211038395)] ),
    ( 'Bax_BclxL_RA' ,
      [Parameter('kbaxbclxl_RAf', 3.057254e-04 / cyto_size  ),
       Parameter('kbaxbclxl_RAr', 0.0036839215247003671),
       Parameter('kbaxbclxl_RAc', 0.15258159856459949)] ),
    ( 'BID_BAX_BAK_inh' ,
      [[Parameter('kbidbcl2f', 1.154239e-05 / mitoM_size  ), Parameter('kbidbcl2r', 0.0012910108044422075)],
       [Parameter('baxbcl2f', 8.587185e-07 / mitoM_size  ), Parameter('baxbcl2r', 0.0065202845242796782)],
       [Parameter('bidbclxlf', 9.986740e-08 / mitoM_size  ), Parameter('bidbclxlr', 0.00041535543108724181)],
       [Parameter('baxbclxlf', 2.886784e-05 / mitoM_size  ), Parameter('baxbclxlr', 0.00025012488950092962)],
       [Parameter('bakbclxlf', 1.780274e-06 / mitoM_size  ), Parameter('bakbclxlr', 0.0004269715362101447)],
       [Parameter('bidmcl1f', 4.710147e-05 / mitoM_size  ), Parameter('bidmcl1r', 0.00040766244776426321)],
       [Parameter('bakmcl1f', 3.248383e-06 / mitoM_size  ), Parameter('bakmcl1r', 0.00013639842040096307)]] ),
    ( 'BCLs_sens' ,
      [[Parameter('kbadbcl2f', 3.433972e-06 / mitoM_size  ), Parameter('kbadbcl2r', 0.00022704112768461186)],
       [Parameter('kbadbclxlf', 1.059135e-06 / mitoM_size  ), Parameter('kbadbclxlr', 0.00077968725403778337)],
       [Parameter('knoxamcl1f', 6.385763e-06 / mitoM_size  ), Parameter('knoxamcl1r', 0.00012652961118041968)]] ),
    ( 'BAX_CYTC' ,
      [[Parameter('kbaxcytocMCf', 3.411127e-05 / mito_size  ),
        Parameter('kbaxcytocMCr', 0.00064859877429990155),
        Parameter('kbaxcytocMCc', 1.2592069541452351)]] ),
    ( 'BAX_SMAC' ,
      [[Parameter('kbaxsmacCAf', 2.326482e-06 / mito_size  ),
        Parameter('kbaxsmacCAr', 5.5807132501168955e-05),
        Parameter('kbaxsmacCAc', 3.1955587174649978)]] ),
    ( 'BAK_CYTC' ,
      [[Parameter('kbakcytocMCf', 4.845554e-05 / mito_size  ),
        Parameter('kbakcytocMCr', 3.7167944717555268e-05),
        Parameter('kbakcytocMCc', 11.46834505509575)]] ),
    ( 'BAK_SMAC' ,
      [[Parameter('kbaksmacCAf', 2.139723e-05 / mito_size  ),
        Parameter('kbaksmacCAr', 0.0001097322748210234),
        Parameter('kbaksmacCAc', 3.9038699062957178)]] ),
    ( 'CYTOC_ACT' ,
      [Parameter('kcytocCcytoAf', 0.0006877586817236283), Parameter('kcytocCcytoAr', 0.0015503474102571645)] ),
    ( 'SMAC_ACT' ,
      [Parameter('ksmacCsmacAf', 0.041134274415515475), Parameter('ksmacCsmacAr', 0.00035726399778643586)] ),

    # PORE TO PARP Parameters
    # -----------------------
    #
    ( 'APAF_CYTC' ,
      [Parameter('kcytocCapaff', 7.149498e-07 / cyto_size  ),
       Parameter('kcytocCapafr', 0.00089176263305991519),
       Parameter('kcytocCapafc', 1.0223679718831065)] ),
    ( 'APOP_C9:APAF' ,
      [Parameter('kapafc9f', 5.960912e-08 / cyto_size  ), Parameter('kapafc9r', 0.00097329314613755786)] ),
    ( 'APOP_C3' ,
      [Parameter('kapopc3f', 6.178960e-09 / cyto_size  ),
       Parameter('kapopc3r', 0.00041079476241540991),
       Parameter('kapopc3c', 1.3489060960433075)] ),
    ( 'APOP_XIAP' ,
      [Parameter('kapopxiapf', 1.584829e-06 / cyto_size  ), Parameter('kapopxiapr', 0.00096233680209892325)] ),
    ( 'SMAC_XIAP' ,
      [Parameter('ksmacxiapf', 7.254456e-06 / cyto_size  ), Parameter('ksmacxiapr', 0.0012839862735030845)] ),
    ( 'C3_C8' ,
      [Parameter('kc8c3f', 1.221922e-07 / cyto_size  ),
       Parameter('kc8c3r', 0.00081533326999884848),
       Parameter('kc8c3c', 0.60873311638624716)] ),
    ( 'C6_C3' ,
      [Parameter('kc3c6f', 1.447397e-06 / cyto_size  ),
       Parameter('kc3c6r', 0.00091981590548697255),
       Parameter('kc3c6c', 0.76370852790888122)] ),
    ( 'C3_XIAP' ,
      [Parameter('kxiapc3f', 1.540781e-06 / cyto_size  ),
       Parameter('kxiapc3r', 0.00044874286581983956),
       Parameter('kxiapc3c', 0.11287003069265149)] ),
    ( 'PARP_C3' ,
      [Parameter('kc3parpf', 1.208152e-06 / cyto_size  ),
       Parameter('kc3parpr', 0.0096256305667117396),
       Parameter('kc3parpc', 1.4375617436356576)] ),
    ( 'C8_C6' ,
      [Parameter('kc6c8f', 5.510615e-08 / cyto_size  ),
       Parameter('kc6c8r', 0.001136332879789306),
       Parameter('kc6c8c', 1.0603741861187763)] ),
    ( 'INIT_AMTS' ,
      [ Parameter('Trail_0'  ,  3000), # 3000 Ligand correspond to 50 ng/ml SuperKiller TRAIL
        Parameter('DR4_0'    ,   200), # 200 death receptor 4
        Parameter('DR5_0'    ,   200), # 200 death receptor 5
        Parameter('Fadd_0'   , 1.0e3), # Fadd
        Parameter('flip_0'   , 1.0e2), # Flip
        Parameter('C8_0'     , 2.0e4), # procaspase-8 
        Parameter('BAR_0'    , 1.0e3), # Bifunctional apoptosis regulator
        Parameter('Bid_0'    , 4.0e4), # Bid
        Parameter('Bax_0'    , 0.8e5), # Bax
        Parameter('Bak_0'    , 0.2e5), # Bak
        Parameter('Bcl2_0'   , 2.0e4), # cytosolic Bcl2
        Parameter('BclxL_0'  , 2.0e4), # cytosolic BclxL
        Parameter('Mcl1_0'   , 2.0e4), # mitochondrial Mcl1  
        Parameter('Bad_0'    , 1.0e3), # Bad
        Parameter('NOXA_0'   , 1.0e3), # NOXA
        Parameter('CytoC_0'  , 5.0e5), # cytochrome c
        Parameter('Smac_0'   , 1.0e5), # Smac    
        Parameter('Apaf_0'   , 1.0e5), # Apaf-1
        Parameter('C3_0'     , 1.0e4), # procaspase-3 (pro-C3)
        Parameter('C6_0'     , 1.0e4), # procaspase-6 (pro-C6)  
        Parameter('C9_0'     , 1.0e5), # procaspase-9 (pro-C9)
        Parameter('XIAP_0'   , 1.0e5), # X-linked inhibitor of apoptosis protein  
        Parameter('PARP_0'   , 1.0e6), # C3* substrate
        ])
    ])


