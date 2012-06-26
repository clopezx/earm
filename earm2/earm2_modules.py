from pysb import *
from pysb.util import alias_model_components
from earm2.macros import *

transloc_rates = [0.01, 0.01]
bcl2_rates = [1.428571e-05, 1e-3] # #1.0e-6/v
bid_effector_rates = [1e-7, 1e-3, 1] # Generalize to catalysis rates?

def earm2_pores():
    alias_model_components()

    # Rates
    pore_max_size = 4
    pore_rates = [[2.040816e-04,  # 1.0e-6/v**2
                   1e-3]] * (pore_max_size - 1)
    pore_transport_rates = [2.857143e-5, 1e-3, 10] # 2e-6 / v?

    # Pore formation by effectors
    # ------------------------------------
    assemble_pore_sequential(Bax(bf=None, state='A'), pore_max_size, pore_rates)
    assemble_pore_sequential(Bak(bf=None, state='A'), pore_max_size, pore_rates)

    # CytC, Smac release
    # ----------------------
    # ringp_transport(Subunit, Source, Dest, min_size, max_size, rates):
    pore_transport(Bax, 4, 4, CytoC(state='M'), CytoC(state='C'),
                   pore_transport_rates)
    pore_transport(Bax, 4, 4, Smac(state='M'), Smac(state='C'),
                   pore_transport_rates)
    pore_transport(Bak, 4, 4, CytoC(state='M'), CytoC(state='C'),
                   pore_transport_rates)
    pore_transport(Bak, 4, 4, Smac(state='M'), Smac(state='C'),
                   pore_transport_rates)

    # --------------------------------------
    # CytC and Smac activation after release
    # --------------------------------------
    equilibrate(Smac(bf=None, state='C'), Smac(bf=None, state='A'),
                          transloc_rates)

    equilibrate(CytoC(bf=None, state='C'), CytoC(bf=None, state='A'),
                          transloc_rates)

def embedded():
    alias_model_components()
    # Args to describe Bax or Bak in the active state but not
    # yet oligomerized

    # tBID to MOMP 
    # ======================
    # tBid and Bad spontaneously insert into the mito membrane
    # --------------------------------------------------------
    # TODO: Bid translocation rates?
    equilibrate(Bid(bf=None, state='T'), Bid(bf=None, state='M'),
                transloc_rates)
    equilibrate(Bad(bf=None, state='C'), Bad(bf=None, state='M'),
                transloc_rates)

    # tBid in the mitochondria activates Bax(C) and Bak(M)
    #-----------------------------------------------------
    catalyze(Bid(state='M'), Bax(state='C'), Bax(state='A'), bid_effector_rates)
    catalyze(Bid(state='M'), Bak(state='M'), Bak(state='A'), bid_effector_rates)

    # Autoactivation, Bax and Bak activate their own kind, but only when
    # free (i.e. not part of a pore complex)
    # ------------------------------------------------------------------
    effector_auto_rates = [1e-7, 1e-3, 1]
    active_monomer = {'s1': None, 's2': None, 'state': 'A'}
    catalyze(Bax(**active_monomer), Bax(state='C'), Bax(state='A'),
             effector_auto_rates)
    catalyze(Bak(**active_monomer), Bak(state='M'), Bak(state='A'),
             effector_auto_rates)

    # ------------------------------------
    # MOMP Inhibition
    # -----------------------------------------------------
    # tBid and free Bax recruit Bcl-xL(C)
    # ------------------------------------------------------------------------
    bclxl_recruitment_rates = [2.040816e-04, 1e-3, 1]
    catalyze(Bid(state='M'), BclxL(state='C'), BclxL(state='M'),
             bclxl_recruitment_rates)
    catalyze(Bax(**active_monomer), BclxL(state='C'), BclxL(state='M'),
             bclxl_recruitment_rates)

    # Bcl2 inhibitors of Bax, Bak, and Bid
    # ------------------------------------
    bind_table([[                              Bcl2,  BclxL(state='M'),        Mcl1],
                [Bid(state='M'),         bcl2_rates,        bcl2_rates,  bcl2_rates],
                [Bax(**active_monomer),  bcl2_rates,        bcl2_rates,        None],
                [Bak(**active_monomer),        None,        bcl2_rates,  bcl2_rates]])
                      
    # Bcl2 sensitizers 
    # ---------------------------------------------------------------
    bind_table([[                       Bcl2,  BclxL(state='M'),        Mcl1],
                [Bad(state='M'),  bcl2_rates,        bcl2_rates,        None],
                [NOXA,                  None,              None,  bcl2_rates]], 


def indirect():
    alias_model_components()
    # tBID to MOMP 
    # ======================
    # BclxL, Bid, Bax migration to mitochondria
    # ----------------------------------------
    # TODO: Bid doesn't migrate? Should BclXL migrate?
    free_Bax = Bax(bf=None, s1=None, s2=None)
    equilibrate(free_Bax(state='C'), free_Bax(state='A'),
                transloc_rates)
    equilibrate(BclxL(bf=None, state='C'), BclxL(bf=None, state='M'),
                transloc_rates)

    # ------------------------------------
    # MOMP Inhibition
    # ------------------------------------
    # Bcl2 inhibitors of Bax, Bak, and Bid
    # ------------------------------------
    #NOTE: indirect not clear about state of Bcl-xL
    bind_table([[                                   BclxL(state='M'),        Mcl1],
                [Bid(state='T'),                          bcl2_rates,  bcl2_rates],
                [Bax(s1=None, s2=None, state='A'),        bcl2_rates,       False],
                [Bak(s1=None, s2=None, state='A'),        bcl2_rates,  bcl2_rates]],

    # Bcl2 sensitizers 
    # ---------------------------------------------------------------
    #Note: according to Andrews, these should all be at mitochondria
    bind_table([[       BclxL(state='M'),        Mcl1],
                [Bad,         bcl2_rates,       False],
                [NOXA,             False,  bcl2_rates]])

         
def direct():
    alias_model_components()
    # tBID to MOMP 
    # ======================
    # Bid, Bax, BclxL migration to mitochondria
    # ----------------------------------------
    # NOTE non-default translocation rates for tBid, suggesting equilibrium
    # strongly favoring membrane binding
    equilibrate(Bid(bf=None, state='T'), Bid(bf=None, state='M'), [1e-1, 1e-3])

    free_Bax = Bax(bf=None, s1=None, s2=None) # Alias for readability
    equilibrate(free_Bax(state='C'), free_Bax(state='M'), transloc_rates)

    equilibrate(BclxL(bf=None, state='C'), BclxL(bf=None, state='M'),
                transloc_rates)

    # Mitochondrial or Cytosolic tBid activates Bax/Bak
    # Bax/Bak form pores
    # ------------------------------------
    catalyze(Bid(state='M'), Bax(state='M'), Bax(bf=None, state='A'),
             bid_effector_rates)
    catalyze(Bid(state='M'), Bak(state='M'), Bak(bf=None, state='A'),
             bid_effector_rates)

    # ------------------------------------
    # MOMP Inhibition
    # ------------------------------------
    # Bcl2 inhibitors of Bax, Bak, and Bid
    # ------------------------------------
    # TODO: Bax/Bak not inhibited in direct model
    bind_table([[                       Bcl2, BclxL(state='M'),        Mcl1],
                [Bid(state='M'),  bcl2_rates,       bcl2_rates,  bcl2_rates]])

    # Bcl2 sensitizers
    # ---------------------------------------------------------------
    bind_table([[             Bcl2,  BclxL(state='M'),        Mcl1],
                [Bad,   bcl2_rates,        bcl2_rates,        None],
                [NOXA,        None,              None,  bcl2_rates]])



