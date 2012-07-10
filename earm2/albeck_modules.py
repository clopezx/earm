"""EARM 1.0 MODULES

#
# These segments are adapted from:
# Albeck JG, Burke JM, Spencer SL, Lauffenburger DA, Sorger PK, 2008
# Modeling a Snap-Action, Variable-Delay Switch Controlling Extrinsic
# Cell Death. PLoS Biol 6(12): e299. doi:10.1371/journal.pbio.0060299
#
# http://www.plosbiology.org/article/info:doi/10.1371/journal.pbio.0060299
"""

from pysb import *
from pysb.util import alias_model_components
from earm2.macros import *
from pysb.macros import equilibrate

# Default forward, reverse, and catalytic rates
KF = 1e-6
KR = 1e-3
KC = 1

#
# RECEPTOR TO BID SEGMENT
def rec_to_bid():
    """Defines the interactions from the ligand insult (e.g. TRAIL) to Bid
    activation as per EARM 1.0.

    This function depends specifically
    on the parameters and monomers of earm_1_0 to work. This function
    uses L, R, DISC, flip, C8, BAR, and Bid monomers and their
    associated parameters to generate the rules that describe Ligand
    to Receptor binding, DISC formation, Caspase8 activation and
    inhibition by flip and BAR as specified in EARM1.0.
    """

    # Needed to recognize the monomer and parameter names in the present scope
    alias_model_components()

    # =====================
    # tBID Activation Rules
    # ---------------------
    #        L + R <--> L:R --> DISC
    #        pC8 + DISC <--> DISC:pC8 --> C8 + DISC
    #        Bid + C8 <--> Bid:C8 --> tBid + C8
    # ---------------------
    two_step_conv(L(), R(), DISC(bf=None ), [4e-7, KR, 1e-5])
    catalyze(DISC(), C8(state='pro'), C8(state='A'), [KF, KR, KC])
    catalyze(C8(state='A'), Bid(state='U'), Bid(state='T'), [KF, KR, KC])
    # ---------------------
    # Inhibition Rules
    # ---------------------
    #        flip + DISC <-->  flip:DISC  
    #        C8 + BAR <--> BAR:C8 CSPS
    # ---------------------
    bind(DISC(), flip(), [KF, KR])
    bind(BAR(), C8(state='A'), [KF, KR])

def pore_to_parp():
    """ This module defines what happens after the pore is activated and CytC
    and Smac are released.

    This is a very specific function which depends on specifically
    on the parameters and monomers of earm_1_0 to work. This function
    uses, CytoC, Smac, Apaf, Apop, C3, C6, C8, C9, PARP, XIAP
    monomers and their associated parameters to generate the rules
    that describe CytC and Smac export from the mitochondria by the
    active pore activation of Caspase3, loopback through Caspase 6,
    and some inhibitions as specified in EARM1.0.
    """

    alias_model_components()
    # ========================
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
    bind(Smac(state='A'), XIAP(), [7e-6, KR]) 
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
    catalyze(C3(state='A'), PARP(state='U'), PARP(state='C'), [KF, 1e-2, KC])
    catalyze(C3(state='A'), C6(state='pro'), C6(state='A'), [KF, KR, KC])
    catalyze(C6(state='A'), C8(state='pro'), C8(state='A'), [3e-8, KR, KC])

def momp_monomers():
    # == Activators ===================
    # Bid, states: Untruncated, Truncated, truncated and Mitochondrial
    Monomer('Bid', [site_name, 'state'], {'state':['U', 'T', 'M']})
    # == Effectors ====================
    # Bax, states: Cytoplasmic, Mitochondrial, Active
    # sites 's1' and 's2' are used for pore formation
    Monomer('Bax', [site_name, 's1', 's2', 'state'], {'state':['C', 'M', 'A']})
    # == Anti-Apoptotics ==============
    Monomer('Bcl2', [site_name])

    # == Cytochrome C and Smac ========
    Monomer('CytoC', [site_name, 'state'], {'state':['M', 'C', 'A']})
    Monomer('Smac', [site_name, 'state'], {'state':['M', 'C', 'A']})

def Bax_tetramerizes(bax_active_state='A', rate_scaling_factor=1):
    active_unbound = {'state': bax_active_state, 'bf': None}
    active_bax_monomer = Bax(s1=None, s2=None, **active_unbound)
    bax2 =(Bax(s1=1, s2=None, **active_unbound) %
           Bax(s1=None, s2=1, **active_unbound))
    bax4 =(Bax(s1=1, s2=4, **active_unbound) %
           Bax(s1=2, s2=1, **active_unbound) %
           Bax(s1=3, s2=2, **active_unbound) %
           Bax(s1=4, s2=3, **active_unbound))
    Rule('Bax_dimerization', active_bax_monomer + active_bax_monomer <> bax2,
         Parameter('Bax_dimerization_kf', KF*rate_scaling_factor),
         Parameter('Bax_dimerization_kr', KR))
    Rule('Bax_tetramerization', bax2 + bax2 <> bax4,
         Parameter('Bax_tetramerization_kf', 2*KF*rate_scaling_factor),
         Parameter('Bax_tetramerization_kr', 0.5*KR))

def Bcl2_binds_Bax1_Bax2_and_Bax4(bax_active_state='A', rate_scaling_factor=1):
    bind(Bax(state=bax_active_state, s1=None, s2=None), Bcl2,
         [KF*rate_scaling_factor, KR])
    pore_bind(Bax(state=bax_active_state), 's1', 's2', 'bf', 2, Bcl2, 'bf',
         [KF*rate_scaling_factor, KR])
    pore_bind(Bax(state=bax_active_state), 's1', 's2', 'bf', 4, Bcl2, 'bf',
         [KF*rate_scaling_factor, KR])

## MOMP Module Implementations
def albeck_11b(do_pore_transport=False):
    """ TODO: Docstring """
    alias_model_components()

    # Set initial conditions
    Initial(Bid(state='U', bf=None), Parameter('Bid_0', 1e5))
    Initial(Bax(state='C', s1=None, s2=None, bf=None), Parameter('Bax_0', 1e5))
    Initial(Bcl2(bf=None), Parameter('Bcl2_0', 2e4))

    # MOMP Mechanism
    catalyze(Bid(state='T'), Bax(state='C', s1=None, s2=None),
             Bax(state='A', s1=None, s2=None), [1e-7, KR, KC])
    bind(Bid(state='T'), Bcl2, [0, KR])
    bind(Bax(state='A', s1=None, s2=None), Bcl2, [KF, KR])

    # Transport of Smac and Cytochrome C
    if do_pore_transport:
        Initial(Smac(state='M', bf=None), Parameter('Smac_0', 1e6))
        #Initial(CytoC(state='M', bf=None), Parameter('CytoC_0', 1e6))
        catalyze(Bax(state='A'), Smac(state='M'), Smac(state='C'),
            [KF, KR, 10])
        #catalyze(Bax(state='A'), CytoC(state='M'), CytoC(state='C'),
        #    [KF, KR, 10])

def albeck_11c(do_pore_transport=False):
    """ TODO: Docstring """
    alias_model_components()
    Initial(Bid(state='U', bf=None), Parameter('Bid_0', 4e4))
    Initial(Bax(state='C', s1=None, s2=None, bf=None), Parameter('Bax_0', 1e5))
    Initial(Bcl2(bf=None), Parameter('Bcl2_0', 2e4))

    # tBid activates Bax
    catalyze(Bid(state='T'), Bax(state='C', s1=None, s2=None),
             Bax(state='A', s1=None, s2=None), [1e-7, KR, KC])

    # Bax tetramerizes
    Bax_tetramerizes(bax_active_state='A')

    # Bcl2 inhibits Bax, Bax2, and Bax4
    Bcl2_binds_Bax1_Bax2_and_Bax4(bax_active_state='A')

    if do_pore_transport:
        Initial(Smac(state='M', bf=None), Parameter('Smac_0', 1e6))
        #Initial(CytoC(state='M', bf=None), Parameter('CytoC_0', 1e6))
        # NOTE change in KF here from previous model!!!!
        pore_transport(Bax(state='A'), Smac(state='M'), Smac(state='C'),
            [[2*KF, KR, 10]])
        #pore_transport(Bax(state='A'), CytoC(state='M'), CytoC(state='C'),
        #    [[KF, KR, 10]])

# Needs separate mitochondrial compartment--by rate scaling?
def albeck_11d(do_pore_transport=False):
    alias_model_components()
    Initial(Bid(state='U', bf=None), Parameter('Bid_0', 4e4))
    Initial(Bax(state='C', s1=None, s2=None, bf=None), Parameter('Bax_0', 1e5))
    Initial(Bcl2(bf=None), Parameter('Bcl2_0', 2e4))

    v = 0.07
    rate_scaling_factor = 1./v

    # tBid activates Bax in the cytosol
    catalyze(Bid(state='T'), Bax(state='C', s1=None, s2=None),
             Bax(state='A', s1=None, s2=None), [1e-7, KR, KC])

    # Active Bax translocates to the mitochondria
    equilibrate(Bax(state='A', bf=None, s1=None, s2=None),
                Bax(state='M', bf=None, s1=None, s2=None),
                [1e-2, 1e-2])

    # Bax tetramerizes
    Bax_tetramerizes(bax_active_state='M',
                     rate_scaling_factor=rate_scaling_factor)

    # Bcl2 inhibits Bax, Bax2, and Bax4
    Bcl2_binds_Bax1_Bax2_and_Bax4(bax_active_state='M',
                                  rate_scaling_factor=rate_scaling_factor)

    if do_pore_transport:
        Initial(Smac(state='M', bf=None), Parameter('Smac_0', 1e6))
        #Initial(CytoC(state='M', bf=None), Parameter('CytoC_0', 1e6))
        # NOTE change in KF here from previous model!!!!
        pore_transport(Bax(state='M'), Smac(state='M'), Smac(state='C'),
            [[rate_scaling_factor*2*KF, KR, 10]])
        #pore_transport(Bax(state='A'), CytoC(state='M'), CytoC(state='C'),
        #    [[KF, KR, 10]])

def albeck_11e(do_pore_transport=False):
    """ TODO: Docstring """
    # Build off of the previous model
    albeck_11d(do_pore_transport=False)

    # Add the "Mito" species, with states "Inactive" and "Active".
    Monomer('Mito', ['bf', 'state'], {'state': ['I', 'A']})
    alias_model_components()
    Initial(Mito(state='I', bf=None), Parameter('Mito_0', 5e5))

    v = 0.07
    rate_scaling_factor = 1./v

    # Add Binding of Bax4 to Mito
    pore_bind(Bax(state='M'), 's1', 's2', 'bf', 4, Mito(state='I'), 'bf',
         [KF*rate_scaling_factor, KR])
    Rule('Mito_activation',
         MatchOnce(Bax(state='M', bf=5, s1=1, s2=4) %
                   Bax(state='M', bf=None, s1=2, s2=1) %
                   Bax(state='M', bf=None, s1=3, s2=2) %
                   Bax(state='M', bf=None, s1=4, s2=3) %
                   Mito(state='I', bf=5)) >>
                   Mito(state='A', bf=None),
         Parameter('Mito_activation_kc', KC))

    if do_pore_transport:
        Initial(Smac(state='M', bf=None), Parameter('Smac_0', 1e6))
        #Initial(CytoC(state='M', bf=None), Parameter('CytoC_0', 1e6))
        catalyze(Mito(state='A'), Smac(state='M'), Smac(state='C'),
            [rate_scaling_factor*2*KF, KR, 10])
        #pore_transport(Bax(state='A'), CytoC(state='M'), CytoC(state='C'),
        #    [[KF, KR, 10]])

def albeck_11f(do_pore_transport=False):
    albeck_11e(do_pore_transport=do_pore_transport)
    alias_model_components()

    # Set parameter values for cooperative pore formation
    equilibrate_BaxA_to_BaxM_kf.value = 1e-4  # was 1e-2 in 11e
    equilibrate_BaxA_to_BaxM_kr.value = 1e-4  # was 1e-2 in 11e
    Bax_dimerization_kf.value /= 100          # was 1e-6 in 11e
    Bax_tetramerization_kf.value /= 10      # was 1e-6 in 11e


