from pysb import *
from collections import OrderedDict

# EARM 2.0 embedded parameters
# from simulated annealing fits to experimental data using a sum of chi-sq objective functions
#

# Reaction rates
parameter_dict = OrderedDict([
    # Trail to Bid parameters
    # -----------------------
    # rec_to_bid module parameters
    # TT binding to DR monomer- 4,5
   ('TT_DRmono',  [ Parameter('ttdrmf', 4.0e-07),       
                    Parameter('ttdrmr', 1.0e-03)]),
    # TT binding to DR dimers- 4,5
   ('TT_DRdim',   [ Parameter('ttdrdf', 4.0e-07),     
                    Parameter('ttdrdr', 1.0e-03)]),
    # TT binding to DR trimers- 4,5
   ('TT_DRtrim',  [ Parameter('ttdrtf', 4.0e-07),     
                    Parameter('ttdrtr', 1.0e-03)]),
    # DR4 trimerization
   ('DR4_RINGP',  [[Parameter('kdr4dimf',  2.040816e-04), #1.0e-6/v**2
                    Parameter('kdr4dimr',  1.0e-3      )],
                    [Parameter('kdr4trimf', 2.040816e-04),
                    Parameter('kdr4trimr', 1.0e-3      )]]),
    # DR5 trimerization
   ('DR5_RINGP',  [[Parameter('kdr5dimf',  2.040816e-04), #1.0e-6/v**2
                    Parameter('kdr5dimr',  1.0e-3      )],
                   [Parameter('kdr5trimf', 2.040816e-04),
                    Parameter('kdr5trimr', 1.0e-3      )]]),
    # Fadd LDRC complex
   ('Fadd_LDRC',  [ Parameter('faddldrcf', 4.0e-07),       
                    Parameter('faddldrcr', 1.0e-03)]),
    # pC8 FADD binding pC8_fadd_b
   ('pC8_fadd_b', [ Parameter('pc8faddf', 4.0e-07),       
                    Parameter('pc8faddr', 1.0e-03)]),
    # pC8 activation to C8 within the same FADD trimer
   ('pC8_dim_act_s', [ Parameter('pc8dimacts', 1.0e-3)]),
    # pC8 activation to C8 across separate FADD trimers
   ('pC8_dim_act_o', [ Parameter('pc8dimacto', 1.0e-3)]),
    # Bid cleavage by C8
   ('C8_BID',     [ Parameter('kc8bidf', 1.0e-07),
                    Parameter('kc8bidr', 1.0e-03),
                    Parameter('kc8bidc', 1.0    )]),
    # DISC inhibition by FLIP
   ('DISC_FLIP',  [ Parameter('kflipdiscf', 1.0e-06), 
                    Parameter('kflipdiscr', 1.0e-03)]),
    # C8 inhibition by BAR
   ('BAR_C8',     [ Parameter('kbarc8f', 1.0e-06), 
                    Parameter('kbarc8r', 1.0e-03)]),

    # BID TO MOMP parameters
    # ----------------------
    #
                    
    ( 'BID_trans' ,
      [Parameter(name='kbidTbidMf', value=0.0026601962092430266), Parameter(name='kbidTbidMr', value=6.4294317349802732e-05)] ),
    ( 'BAD_trans' ,
      [Parameter(name='kbadCbadMf', value=0.0029419710017770154), Parameter(name='kbadCbadMr', value=0.0027381917924608314)] ),
    ( 'BAX_trans' ,
      [Parameter(name='kbaxCbaxMf', value=0.0027070779785728404), Parameter(name='kbaxCbaxMr', value=0.004353093828158041)] ),
    ( 'BCL2_trans' ,
      [Parameter(name='kbcl2Cbcl2Mf', value=0.0014737610524655136), Parameter(name='kbcl2Cbcl2Mr', value=0.014906933079643386)] ),
    ( 'BCLXL_trans' ,
      [Parameter(name='kbclxlCbclxlMf', value=0.027466171257673082), Parameter(name='kbclxlCbclxlMr', value=0.0084157689975014327)] ),
    ( 'BID_BAX' ,
      [Parameter(name='kbidbaxf', value=3.8956533050487824e-07),
       Parameter(name='kbidbaxr', value=0.00036856169047138091),
       Parameter(name='kbidbaxc', value=1.4664186386875186)] ),
    ( 'BID_BAK' ,
      [Parameter(name='kbidbakf', value=6.1557947201966136e-09),
       Parameter(name='kbidbakr', value=0.0029106994159151453),
       Parameter(name='kbidbakc', value=0.84528236986831273)] ),
    ( 'BAX_BAX' ,
      [Parameter(name='kbaxbaxf', value=3.0468256843466471e-10),
       Parameter(name='kbaxbaxr', value=1.7249510431983885e-05),
       Parameter(name='kbaxbaxc', value=2.1191713445212592)] ),
    ( 'BAK_BAK' ,
      [Parameter(name='kbakbakf', value=3.6635565426345757e-08),
       Parameter(name='kbakbakr', value=0.0013849096067022136),
       Parameter(name='kbakbakc', value=0.23567721868317815)] ),
    ( 'BAX_PORE' ,
      [[Parameter(name='kbaxdimf', value=0.00019873556680652295), Parameter(name='kbaxdimr', value=0.00051552975969350876)],
       [Parameter(name='kbaxtrimf', value=0.00026743128681617405), Parameter(name='kbaxtrimr', value=0.00053921595886761769)],
       [Parameter(name='kbaxtetf', value=5.3331379707057153e-05), Parameter(name='kbaxtetr', value=1.6982589396135685e-05)]] ),
    ( 'BAK_PORE' ,
      [[Parameter(name='kbakdimf', value=0.00044382361842463093), Parameter(name='kbakdimr', value=0.0017427126140823314)],
       [Parameter(name='kbaktrimf', value=0.00073068817590211152), Parameter(name='kbaktrimr', value=0.00041757049236862946)],
       [Parameter(name='kbaktetf', value=1.5857586148331407e-05), Parameter(name='kbaktetr', value=0.0010660677426411207)]] ),
    ( 'Bid_BclxL_RA' ,
      [Parameter(name='kbidbclxl_RAf', value=0.00069342164498488337),
       Parameter(name='kbidbclxl_RAr', value=0.00083046693481875049),
       Parameter(name='kbidbclxl_RAc', value=0.69567061211038395)] ),
    ( 'Bax_BclxL_RA' ,
      [Parameter(name='kbaxbclxl_RAf', value=0.00030572537752843616),
       Parameter(name='kbaxbclxl_RAr', value=0.0036839215247003671),
       Parameter(name='kbaxbclxl_RAc', value=0.15258159856459949)] ),
    ( 'BID_BAX_BAK_inh' ,
      [[Parameter(name='kbidbcl2f', value=1.1542394529607694e-05), Parameter(name='kbidbcl2r', value=0.0012910108044422075)],
       [Parameter(name='baxbcl2f', value=8.587184549762698e-07), Parameter(name='baxbcl2r', value=0.0065202845242796782)],
       [Parameter(name='bidbclxlf', value=9.986740023291606e-08), Parameter(name='bidbclxlr', value=0.00041535543108724181)],
       [Parameter(name='baxbclxlf', value=2.8867841570771106e-05), Parameter(name='baxbclxlr', value=0.00025012488950092962)],
       [Parameter(name='bakbclxlf', value=1.7802744769530515e-06), Parameter(name='bakbclxlr', value=0.0004269715362101447)],
       [Parameter(name='bidmcl1f', value=4.7101474631746376e-05), Parameter(name='bidmcl1r', value=0.00040766244776426321)],
       [Parameter(name='bakmcl1f', value=3.2483828371507278e-06), Parameter(name='bakmcl1r', value=0.00013639842040096307)]] ),
    ( 'BCLs_sens' ,
      [[Parameter(name='kbadbcl2f', value=3.4339717491570382e-06), Parameter(name='kbadbcl2r', value=0.00022704112768461186)],
       [Parameter(name='kbadbclxlf', value=1.0591351543095018e-06), Parameter(name='kbadbclxlr', value=0.00077968725403778337)],
       [Parameter(name='knoxamcl1f', value=6.385762500280488e-06), Parameter(name='knoxamcl1r', value=0.00012652961118041968)]] ),
    ( 'BAX_CYTC' ,
      [[Parameter(name='kbaxcytocMCf', value=3.411127052123396e-05),
        Parameter(name='kbaxcytocMCr', value=0.00064859877429990155),
        Parameter(name='kbaxcytocMCc', value=1.2592069541452351)]] ),
    ( 'BAX_SMAC' ,
      [[Parameter(name='kbaxsmacCAf', value=2.326481927758869e-06),
        Parameter(name='kbaxsmacCAr', value=5.5807132501168955e-05),
        Parameter(name='kbaxsmacCAc', value=3.1955587174649978)]] ),
    ( 'BAK_CYTC' ,
      [[Parameter(name='kbakcytocMCf', value=4.8455542806722519e-05),
        Parameter(name='kbakcytocMCr', value=3.7167944717555268e-05),
        Parameter(name='kbakcytocMCc', value=11.46834505509575)]] ),
    ( 'BAK_SMAC' ,
      [[Parameter(name='kbaksmacCAf', value=2.1397233843771375e-05),
        Parameter(name='kbaksmacCAr', value=0.0001097322748210234),
        Parameter(name='kbaksmacCAc', value=3.9038699062957178)]] ),
    ( 'CYTOC_ACT' ,
      [Parameter(name='kcytocCcytoAf', value=0.0006877586817236283), Parameter(name='kcytocCcytoAr', value=0.0015503474102571645)] ),
    ( 'SMAC_ACT' ,
      [Parameter(name='ksmacCsmacAf', value=0.041134274415515475), Parameter(name='ksmacCsmacAr', value=0.00035726399778643586)] ),

    # PORE TO PARP Parameters
    # -----------------------
    #
    ( 'APAF_CYTC' ,
      [Parameter(name='kcytocCapaff', value=7.1494980862371046e-07),
       Parameter(name='kcytocCapafr', value=0.00089176263305991519),
       Parameter(name='kcytocCapafc', value=1.0223679718831065)] ),
    ( 'APOP_C9:APAF' ,
      [Parameter(name='kapafc9f', value=5.9609124208358764e-08), Parameter(name='kapafc9r', value=0.00097329314613755786)] ),
    ( 'APOP_C3' ,
      [Parameter(name='kapopc3f', value=6.1789602013059776e-09),
       Parameter(name='kapopc3r', value=0.00041079476241540991),
       Parameter(name='kapopc3c', value=1.3489060960433075)] ),
    ( 'APOP_XIAP' ,
      [Parameter(name='kapopxiapf', value=1.5848294920276248e-06), Parameter(name='kapopxiapr', value=0.00096233680209892325)] ),
    ( 'SMAC_XIAP' ,
      [Parameter(name='ksmacxiapf', value=7.2544562746768184e-06), Parameter(name='ksmacxiapr', value=0.0012839862735030845)] ),
    ( 'C3_C8' ,
      [Parameter(name='kc8c3f', value=1.2219219982178813e-07),
       Parameter(name='kc8c3r', value=0.00081533326999884848),
       Parameter(name='kc8c3c', value=0.60873311638624716)] ),
    ( 'C6_C3' ,
      [Parameter(name='kc3c6f', value=1.4473973540999643e-06),
       Parameter(name='kc3c6r', value=0.00091981590548697255),
       Parameter(name='kc3c6c', value=0.76370852790888122)] ),
    ( 'C3_XIAP' ,
      [Parameter(name='kxiapc3f', value=1.5407812547711927e-06),
       Parameter(name='kxiapc3r', value=0.00044874286581983956),
       Parameter(name='kxiapc3c', value=0.11287003069265149)] ),
    ( 'PARP_C3' ,
      [Parameter(name='kc3parpf', value=1.2081523254583362e-06),
       Parameter(name='kc3parpr', value=0.0096256305667117396),
       Parameter(name='kc3parpc', value=1.4375617436356576)] ),
    ( 'C8_C6' ,
      [Parameter(name='kc6c8f', value=5.5106145496793361e-08),
       Parameter(name='kc6c8r', value=0.001136332879789306),
       Parameter(name='kc6c8c', value=1.0603741861187763)] ),
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


