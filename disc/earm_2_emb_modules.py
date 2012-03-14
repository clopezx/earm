# EARM 1.0 MODULES
# Notice the alias_model_components function is needed to recognize the monomer and 
# parameter names in the present scope
#
# rec_to_bid: This module defines the interactions from the ligand insult (e.g. TRAIL)
#             to Bid activation as per EARM 1.0
#
# fxrcsp_to_parp: This module defines what happens after the pore is activated and 
#                 CytC and Smac are released
#
# These segments are adapted from:
# Albeck JG, Burke JM, Spencer SL, Lauffenburger DA, Sorger PK, 2008
# Modeling a Snap-Action, Variable-Delay Switch Controlling Extrinsic
# Cell Death. PLoS Biol 6(12): e299. doi:10.1371/journal.pbio.0060299
#
# http://www.plosbiology.org/article/info:doi/10.1371/journal.pbio.0060299
#
#
from pysbhelperfuncs import *

# get model components accessible in this scope
alias_model_components()

# RECEPTOR TO BID SEGMENT
def bid_to_momp(model, kd):
    """ This is a very specific function which depends on specifically
    on the parameters and monomers of earm_2_0_embedded to work. This function
    uses tBid, Bax, Bak, some inhibitors, and some sensitizers to model the
    Bcl-2 protein family regulation leading to MOMP.
    """
    # tBID to MOMP 
    # ======================
    # tBid and Bad spontaneously insert into the mito membrane
    # --------------------------------------------------------
    Rule('Bid_to_mem', Bid(bf = None, state = 'T') ** cy <> Bid(bf = None, state = 'M') ** mitM, kd['BID_trans'][0],kd['BID_trans'][1])
    Rule('Bad_to_mem', Bad(bf = None, state = 'C') ** cy <> Bad(bf = None, state = 'M') ** mitM, kd['BAD_trans'][0],kd['BAD_trans'][1])
    
    # tBid in the mitochondria activates Bax(C) and Bak(M)
    #-----------------------------------------------------
    two_step_mod(Bid(state = 'M') ** mitM, Bax(state='C') ** cy, Bax(bf = None, state = 'A') ** mitM, kd['BID_BAX'])
    two_step_mod(Bid(state = 'M') ** mitM, Bak(state='M'), Bak(bf = None, state = 'A'), kd['BID_BAK'])
    
    # Autoactivation, Bax and Bak activate their own kind, but only when
    # free (i.e. not part of a pore complex)
    # ------------------------------------------------------------------
    two_step_mod(Bax(state = 'A', s1=None, s2=None) ** mitM, Bax(state='C') ** cy, Bax(bf = None, state = 'A') ** mitM, kd['BAX_BAX'])
    two_step_mod(Bak(state = 'A', s1=None, s2=None), Bak(state='M'), Bak(bf = None, state = 'A'), kd['BAK_BAK'])
    
    # pore_assembly(Subunit, size, rates):
    # ------------------------------------
    ringp_assembly(Bax(bf=None, state='A') ** mitM, 4, kd['BAX_PORE'])
    ringp_assembly(Bak(bf=None, state='A'), 4, kd['BAK_PORE'])
    
    # ------------------------------------
    # MOMP Inhibition
    # -----------------------------------------------------
    # tBid and free Bax recruit Bcl-xL(C) and autoinhibit themselves
    #        Bid/Bax + BclxL <> Bid/Bax:BclxL(C) >> Bid/Bax:BclxL(M)
    #        Notice the product of this feeds into the product of the inh rxns
    # ------------------------------------------------------------------------
    two_step_conv(Bid(state = 'M') ** mitM,           BclxL(state='C') ** cy, Bid(bf = 1,             state = 'M')%BclxL(bf = 1, state='M') ** mitM, kd['Bid_BclxL_RA'])
    two_step_conv(Bax(state = 'A', s1=None, s2=None) ** mitM, BclxL(state='C') ** cy, Bax(bf = 1, s1=None, s2=None, state = 'A')%BclxL(bf = 1, state='M') ** mitM, kd['Bax_BclxL_RA'])
    
    
    # Bcl2 inhibitors of Bax, Bak, and Bid
    # a set of simple bind reactions:
    #        Inh + Act <--> Inh:Act
    # ------------------------------------
    simple_bind_table([[                                            Bcl2,         BclxL,  Mcl1],
                       [                                              {}, {'state':'M'},    {}],
                       [Bid, {'state':'M'},                        True,          True,  True],
                       [Bax, {'s1':None, 's2':None, 'state':'A'},  True,          True, False],
                       [Bak, {'s1':None, 's2':None, 'state':'A'}, False,          True,  True]],
                      kd['BID_BAX_BAK_inh'], model)
    
    # Sensitizers
    # Bcl2 sensitizers bind through a simple bind resction: 
    #        Inh + Act <--> Inh:Act
    # This goes through the list in row-major order (as it should be)
    # ---------------------------------------------------------------
    simple_bind_table([[                      Bcl2,         BclxL,  Mcl1],
                        [                        {}, {'state':'M'},    {}],
                        [Bad,  {'state':'M'},  True,          True, False],
                        [NOXA, {},            False,         False,  True]], 
                        kd['BCLs_sens'], model)
    
    # FXR CASPASES CLEAVE PARP 
    # ========================
    # CytC, Smac release
    # ----------------------
    #        AMito + mCytoC <-->  AMito:mCytoC --> AMito + ACytoC  
    #        ACytoC <-->  cCytoC
    #        AMito + mSmac <-->  AMito:mSmac --> AMito + ASmac  
    #        ASmac <-->  cSmac
    # ----------------------
    # ringp_transport(Subunit, Source, Dest, min_size, max_size, rates):
    ringp_transport(Bax(bf=None), CytoC(state='M'), CytoC(state='C'), 4, 4, kd['BAX_CYTC']) 
    ringp_transport(Bax(bf=None),  Smac(state='M'),  Smac(state='C'), 4, 4, kd['BAX_SMAC']) 
    ringp_transport(Bak(bf=None), CytoC(state='M'), CytoC(state='C'), 4, 4, kd['BAK_CYTC'])
    ringp_transport(Bak(bf=None),  Smac(state='M'),  Smac(state='C'), 4, 4, kd['BAK_SMAC'])
    # --------------------------------------
    # CytC and Smac activation after release
    # --------------------------------------
    Rule('act_cSmac',  Smac(bf=None, state='C') <> Smac(bf=None, state='A'), kd['SMAC_ACT'][0], kd['SMAC_ACT'][1])
    Rule('act_cCytoC', CytoC(bf=None, state='C') <> CytoC(bf=None, state='A'), kd['CYTOC_ACT'][0], kd['CYTOC_ACT'][1])

def pore_to_parp(model, kd):
    """ This is a very specific function which depends on specifically
    on the parameters and monomers of earm_1_0 to work. This function
    uses, CytoC, Smac, Apaf, Apop, C3, C6, C8, C9, PARP, XIAP
    monomers and their associated parameters to generate the rules
    that describe CytC and Smac export from the mitochondria by the
    active pore activation of Caspase3, loopback through Caspase 6,
    and some inhibitions as specified in EARM1.0.
    """
    # ---------------------------
    # Apoptosome formation
    # ---------------------------
    #        Apaf + cCytoC <-->  Apaf:cCytoC --> aApaf + cCytoC
    #        aApaf + pC9 <-->  Apop
    #        Apop + pC3 <-->  Apop:pC3 --> Apop + C3
    # ---------------------------
    two_step_mod(CytoC(state='A'), Apaf(state='I'), Apaf(bf = None, state = 'A'), kd['APAF_CYTC'])
    one_step_conv(Apaf(state='A'), C9(), Apop(bf=None), kd['APOP_C9:APAF'])
    two_step_mod(Apop(), C3(state='pro'), C3(bf = None, state='A'), kd['APOP_C3'])
    # -----------------------------
    # Apoptosome related inhibitors
    # -----------------------------
    #        Apop + XIAP <-->  Apop:XIAP  
    #        cSmac + XIAP <-->  cSmac:XIAP  
    simple_bind(Apop(), XIAP(), kd['APOP_XIAP'])
    simple_bind(Smac(state='A'), XIAP(), kd['SMAC_XIAP'])
    # ---------------------------
    # Caspase reactions (effectors, inhibitors, and loopback initiators)
    # ---------------------------
    #        pC3 + C8 <--> pC3:C8 --> C3 + C8 CSPS
    #        pC6 + C3 <--> pC6:C3 --> C6 + C3 CSPS
    #        XIAP + C3 <--> XIAP:C3 --> XIAP + C3_U CSPS
    #        PARP + C3 <--> PARP:C3 --> CPARP + C3 CSPS
    #        pC8 + C6 <--> pC8:C6 --> C8 + C6 CSPS
    # ---------------------------
    two_step_mod(C8(state='A'), C3(state='pro'), C3(bf = None, state='A'), kd['C3_C8'])
    two_step_mod(C3(state='A'), C6(state='pro'), C6(bf = None, state='A'), kd['C6_C3'])
    two_step_mod(XIAP(), C3(state = 'A'), C3(bf = None, state = 'ub'), kd['C3_XIAP'])
    two_step_mod(C3(state = 'A'), PARP(state='U'), PARP(bf = None, state='C'), kd['PARP_C3'])
    two_step_mod(C6(state='A'), C8(state='pro'), C8(bf = None, state = 'A'), kd['C8_C6'])
