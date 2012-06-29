"""TODO: Write docstring for earm2_modules

More stuff here.
"""

from pysb import *
from pysb.util import alias_model_components
from earm2.macros import *
from pysb.macros import equilibrate

transloc_rates =     [        1e-2, 1e-2]
bcl2_rates =         [1.428571e-05, 1e-3]    # #1.0e-6/v
bid_effector_rates = [        1e-7, 1e-3, 1] # Generalize to catalysis rates?
# Site-value arguments to describe Bax or Bak in the active state but not
# yet oligomerized
active_monomer = {'s1': None, 's2': None, 'state': 'A'}

def earm2_pore_formation():
    """ Pore formation and transport process used by all modules.
    """

    alias_model_components()

    # Rates
    pore_max_size = 4
    pore_rates = [[2.040816e-04,  # 1.0e-6/v**2
                   1e-3]] * (pore_max_size - 1)
    pore_transport_rates = [[2.857143e-5, 1e-3, 10]] # 2e-6 / v?

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

def translocate_tBid_Bax_BclxL():
    equilibrate(Bid(bf=None, state='T'), Bid(bf=None, state='M'), [1e-1, 1e-3])
    # equilibrate(Bid(bf=None, state='T'), Bid(bf=None, state='M'), transloc_rates)
    # Previous indirect model had more membrane-favorable rates for tBid

    free_Bax = Bax(bf=None, s1=None, s2=None) # Alias for readability
    equilibrate(free_Bax(state='C'), free_Bax(state='M'),
                transloc_rates)

    equilibrate(BclxL(bf=None, state='C'), BclxL(bf=None, state='M'),
                transloc_rates)

def tBid_activates_Bax_and_Bak():
    catalyze(Bid(state='M'), Bax(state='M'), Bax(state='A'), bid_effector_rates)
    catalyze(Bid(state='M'), Bak(state='M'), Bak(state='A'), bid_effector_rates)

def sensitizers_bind_anti_apoptotics():
    bind_table([[                       Bcl2,  BclxL(state='M'),  Mcl1(state='M')],
                [Bad(state='M'),  bcl2_rates,        bcl2_rates,             None],
                [NOXA(state='M'),       None,              None,       bcl2_rates]])

# Embedded Model ===================================================
def embedded():
    """ Direct and indirect modes of action, occurring at the membrane.
    """
    alias_model_components()

    translocate_tBid_Bax_BclxL()

    tBid_activates_Bax_and_Bak()

    # Autoactivation, Bax and Bak activate their own kind, but only when
    # free (i.e. not part of a pore complex)
    effector_auto_rates = [1e-7, 1e-3, 1]
    catalyze(Bax(active_monomer), Bax(state='C'), Bax(state='A'),
             effector_auto_rates)
    catalyze(Bak(active_monomer), Bak(state='M'), Bak(state='A'),
             effector_auto_rates)

    # tBid and free Bax recruit Bcl-xL(C)
    bclxl_recruitment_rates = [2.040816e-04, 1e-3, 1]
    catalyze(Bid(state='M'), BclxL(state='C'), BclxL(state='M'),
             bclxl_recruitment_rates)
    catalyze(Bax(active_monomer), BclxL(state='C'), BclxL(state='M'),
             bclxl_recruitment_rates)

    # Role of activator tBid is to bind all anti-apoptotics
    # NOTE: Doug Green's MODE 1 inhibition
    bind_table([[                            Bcl2,  BclxL(state='M'),  Mcl1(state='M')],
                [Bid(state='M'),       bcl2_rates,        bcl2_rates,       bcl2_rates]])

    # Anti-apoptotics bind activated effectors
    # NOTE: Doug Green's MODE 2 inhibition
    bind_table([[                            Bcl2,  BclxL(state='M'),        Mcl1],
                [Bax(active_monomer),  bcl2_rates,        bcl2_rates,        None],
                [Bak(active_monomer),        None,        bcl2_rates,  bcl2_rates]])

    sensitizers_bind_anti_apoptotics()

    # Bax and Bak form pores by sequential addition
    earm2_pore_formation()

# Indirect Model ===================================================
def indirect():
    """Bax and Bak spontaneously form pores without activation.
       The "activator" tBid binds all of the anti-apoptotics.
    """

    alias_model_components()

    translocate_tBid_Bax_BclxL()

    # Bax and Bak spontaneously become activated
    free_Bax = Bax(bf=None, s1=None, s2=None) # Alias
    free_Bak = Bak(bf=None, s1=None, s2=None) # Alias
    equilibrate(free_Bax(state='M'), free_Bax(state='A'),
                transloc_rates)
    equilibrate(free_Bak(state='M'), free_Bak(state='A'),
                transloc_rates)

    # Role of activator tBid is to bind all anti-apoptotics
    bind_table([[                            Bcl2,  BclxL(state='M'),  Mcl1(state='M')],
                [Bid(state='M'),       bcl2_rates,        bcl2_rates,       bcl2_rates]])

    # Anti-apoptotics bind "active" Bax and Bak to prevent pore formation
    bind_table([[                            Bcl2,  BclxL(state='M'),  Mcl1(state='M')],
                [Bax(active_monomer),  bcl2_rates,        bcl2_rates,             None],
                [Bak(active_monomer),        None,        bcl2_rates,       bcl2_rates]])

    sensitizers_bind_anti_apoptotics()

    # Bax and Bak form pores by sequential addition
    earm2_pore_formation()

# Direct Model ===================================================
def direct():
    """Anti-apoptotics prevent BH3-onlies from activating Bax and Bak.

    Bax and Bak require activation to be able to form pores.
    The anti-apoptotics don't inhibit activated Bax and Bak; their only role
    is to bind BH3-onlies.
    """

    alias_model_components()

    translocate_tBid_Bax_BclxL()

    tBid_activates_Bax_and_Bak()

    # Anti-apoptotics bind and inhibit the activator tBid
    # Doug Green's "Mode 1" inhibition
    bind_table([[                       Bcl2, BclxL(state='M'),  Mcl1(state='M')],
                [Bid(state='M'),  bcl2_rates,       bcl2_rates,       bcl2_rates]])

    # Non-"activator" BH3-onlies (i.e., "sensitizers") bind the anti-apoptotics
    sensitizers_bind_anti_apoptotics()

    # Bax and Bak form pores by sequential addition
    earm2_pore_formation()
