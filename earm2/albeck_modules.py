from pysb import *
from pysb.util import alias_model_components
from earm2.macros import *

# get model components accessible in this scope
KF = 1e-6
KR = 1e-3
KC = 1

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
# RECEPTOR TO BID SEGMENT
def rec_to_bid():
    """ This function depends specifically
    on the parameters and monomers of earm_1_0 to work. This function
    uses L, R, DISC, flip, C8, BAR, and Bid monomers and their
    associated parameters to generate the rules that describe Ligand
    to Receptor binding, DISC formation, Caspase8 activation and
    inhibition by flip and BAR as specified in EARM1.0.
    """

    alias_model_components()
    # RECEPTOR TO tBID
    # =====================
    # tBID Activation Rules
    # ---------------------
    #        L + R <--> L:R --> DISC
    #        pC8 + DISC <--> DISC:pC8 --> C8 + DISC
    #        Bid + C8 <--> Bid:C8 --> tBid + C8
    # ---------------------
    two_step_conv(L(), R(), DISC(bf=None ), [4e-7, KR, 1e-5])
    catalyze(DISC(), C8(state='pro'), C8(state='A'), [KF, KR, KC])
    catalyze(C8(state='A'), Bid(state='U'), Bid(state='T'), [1e-7, KR, KC])
    # ---------------------
    # Inhibition Rules
    # ---------------------
    #        flip + DISC <-->  flip:DISC  
    #        C8 + BAR <--> BAR:C8 CSPS
    # ---------------------
    bind(DISC(), flip(), [KF, KR])
    bind(BAR(), C8(state='A'), [KF, KR])

def pore_to_parp():
    """ This is a very specific function which depends on specifically
    on the parameters and monomers of earm_1_0 to work. This function
    uses, CytoC, Smac, Apaf, Apop, C3, C6, C8, C9, PARP, XIAP
    monomers and their associated parameters to generate the rules
    that describe CytC and Smac export from the mitochondria by the
    active pore activation of Caspase3, loopback through Caspase 6,
    and some inhibitions as specified in EARM1.0.
    """
    alias_model_components()
    # FXR CASPASES CLEAVE PARP 
    # ========================
    #
    # Apoptosome formation
    # ---------------------------
    #        Apaf + cCytoC <-->  Apaf:cCytoC --> aApaf + cCytoC
    #        aApaf + pC9 <-->  Apop
    #        Apop + pC3 <-->  Apop:pC3 --> Apop + C3
    # ---------------------------
    catalyze(CytoC(state='A'), Apaf(state='I'), Apaf(state='A'), [5e-7, KR, KC])
    one_step_conv(Apaf(state='A'), C9(), Apop(bf=None), [5e-8, KR])
    catalyze(Apop(), C3(state='pro'), C3(bf=None, state='A'), [5e-9, KR, KC])
    # -----------------------------
    # Apoptosome related inhibitors
    # -----------------------------
    #        Apop + XIAP <-->  Apop:XIAP  
    #        cSmac + XIAP <-->  cSmac:XIAP  
    bind(Apop(), XIAP(), [2e-6, KR])
    bind(Smac(state='A'), XIAP(), [7e-6, 1e-3])
    # ---------------------------
    # Caspase reactions (effectors, inhibitors, and loopback initiators)
    # ---------------------------
    #        pC3 + C8 <--> pC3:C8 --> C3 + C8 CSPS
    #        pC6 + C3 <--> pC6:C3 --> C6 + C3 CSPS
    #        XIAP + C3 <--> XIAP:C3 --> XIAP + C3_U CSPS
    #        PARP + C3 <--> PARP:C3 --> CPARP + C3 CSPS
    #        pC8 + C6 <--> pC8:C6 --> C8 + C6 CSPS
    # ---------------------------
    catalyze(C8(state='A'), C3(state='pro'), C3(state='A'), [1e-7, KR, KC])
    catalyze(XIAP(), C3(state='A'), C3(state = 'ub'), [2e-6, KR, 1e-1])
    catalyze(C3(state='A'), PARP(state='U'), PARP(state='C'), [1e-6, 1e-2, KC])
    catalyze(C3(state='A'), C6(state='pro'), C6(state='A'), [KF, KR, KC])
    catalyze(C6(state='A'), C8(state='pro'), C8(state='A'), [3e-8, KR, KC])

## MOMP Module Implementations
def albeck2008_11b():
    alias_model_components()
    catalyze(Bid(state='T'), Bax(state='C'), Bax(state='A'), [KF, KR, KC])
    bind(Bax(state='A'), Bcl2, [KF, KR])
    #catalyze(Bax(state='A'), Smac(loc='c'), Smac(loc='r'), [KF, KR, KC])
    
def albeck2008_11c():
    alias_model_components()
    catalyze(Bid(state='T'), Bax(state='C'), Bax(state='A'), [KF, KR, KC])

    dimerize(Bax(state='A'), Bax2, [KF, KR])
    dimerize(Bax2, Bax4, [KF, KR])

    bind_table([[            Bax,      Bax2,     Bax4],
                [Bcl2,  [KF, KR],  [KF, KR],  [KF, KR]]])

    #catalyze(Bax4, Smac(loc='c'), Smac(loc='r'), [KF, KR, KC])

# Needs separate mitochondrial compartment--by rate scaling?
def albeck2008_11d():
    alias_model_components()
    catalyze(Bid(state='T'), Bax(state='C'), Bax(state='A'), [KF, KR, KC])

    dimerize(Bax(state='A'), Bax2, [KF, KR])
    dimerize(Bax2, Bax4, [KF, KR])

    bind_table([[            Bax,      Bax2,     Bax4],
                [Bcl2,  [KF, KR],  [KF, KR],  [KF, KR]]])

    #catalyze(Bax4, Smac(loc='c'), Smac(loc='r'), [KF, KR, KC])
        
def albeck2008_11e():
    alias_model_components()
    catalyze(Bid(state='T'), Bax(state='C'), Bax(state='A'), [KF, KR, KC])

    dimerize(Bax(state='A'), Bax2, [KF, KR])
    dimerize(Bax2, Bax4, [KF, KR])

    bind_and_convert(Bax4, M, Pore, [KF, KR])

    bind_table([[            Bax,      Bax2,     Bax4],
                [Bcl2,  (KF, KR),  (KF, KR),  (KF, KR)]])

    #catalyze(Pore, Smac(loc='c'), Smac(loc='r'), [KF, KR, KC])

def albeck2008_11f():
    alias_model_components()
    pass



