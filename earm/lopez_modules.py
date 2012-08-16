"""TODO: Write docstring for lopez_modules

More stuff here.
"""

from pysb import *
from pysb.util import alias_model_components
from earm.macros import *
from pysb.macros import equilibrate
import albeck_modules

# Default site name for binding
site_name = 'bf'

transloc_rates =     [        1e-2, 1e-2]
bcl2_rates =         [1.428571e-05, 1e-3]    # 1.0e-6/v
bid_effector_rates = [        1e-7, 1e-3, 1] 

# Site-value arguments to describe Bax or Bak in the active state but not
# yet oligomerized
active_monomer = {'s1': None, 's2': None, 'state': 'A'}

# MONOMER DECLARATION MACROS ================================================
def momp_monomers():
    """Declare the monomers for the Bcl-2 family proteins, Cyto c, and Smac.

    The package variable site_name specifies the name of the site to be used
    for all binding reactions (with the exception of Bax and Bak, which have
    additional sites used for oligomerization).

    The 'state' site denotes various localization and/or activity states of a
    Monomer, with 'C' denoting cytoplasmic localization and 'M' mitochondrial
    localization. Most Bcl-2 proteins have the potential for both cytoplasmic
    and mitochondrial localization, with the exceptions of Bak and Bcl-2,
    which are apparently constitutively mitochondrial.
    """

    # == Activators ===================
    # Bid, states: Untruncated, Truncated, truncated and Mitochondrial
    Monomer('Bid', [site_name, 'state'], {'state':['U', 'T', 'M']})
    # == Effectors ====================
    # Bax, states: Cytoplasmic, Mitochondrial, Active
    # sites 's1' and 's2' are used for pore formation
    Monomer('Bax', [site_name, 's1', 's2', 'state'], {'state':['C', 'M', 'A']})
    # Bak, states: inactive and Mitochondrial, Active (and mitochondrial)
    # sites 's1' and 's2' are used for pore formation
    Monomer('Bak', [site_name, 's1', 's2', 'state'], {'state':['M', 'A']})
    # == Anti-Apoptotics ==============
    Monomer('Bcl2', [site_name])
    Monomer('BclxL', [site_name, 'state'], {'state':['C', 'M']})
    Monomer('Mcl1', [site_name, 'state'], {'state':['M', 'C']})
    # == Sensitizers ==================
    Monomer('Bad', [site_name, 'state'], {'state':['C', 'M']})
    Monomer('Noxa', [site_name, 'state'], {'state': ['C', 'M']})

    # == Cytochrome C and Smac ========
    Monomer('CytoC', [site_name, 'state'], {'state':['M', 'C', 'A']})
    Monomer('Smac', [site_name, 'state'], {'state':['M', 'C', 'A']})

# MOMP SEGMENT ==============================================================

## Macros -------------------------------------------------------------------
def declare_initial_conditions():
    """Declare initial conditions for Bcl-2 family proteins, Cyto c, and Smac.
    """
    Parameter('Bid_0'   , 4.0e4) # Bid
    Parameter('BclxL_0' , 2.0e4) # cytosolic BclxL
    Parameter('Mcl1_0'  , 2.0e4) # Mitochondrial Mcl1
    Parameter('Bcl2_0'  , 2.0e4) # Mitochondrial Bcl2
    Parameter('Bad_0'   , 1.0e3) # Bad
    Parameter('Noxa_0'  , 1.0e3) # Noxa
    Parameter('CytoC_0' , 5.0e5) # cytochrome c
    Parameter('Smac_0'  , 1.0e5) # Smac
    Parameter('Bax_0'   , 0.8e5) # Bax
    Parameter('Bak_0'   , 0.2e5) # Bak

    alias_model_components()

    Initial(Bid(bf=None, state='U'), Bid_0)
    Initial(Bad(bf=None, state='C'), Bad_0)
    Initial(Bax(bf=None, s1=None, s2=None, state='C'), Bax_0)
    Initial(Bak(bf=None, s1=None, s2=None, state='M'), Bak_0)
    Initial(Bcl2(bf=None), Bcl2_0)
    Initial(BclxL (bf=None, state='C'), BclxL_0)
    Initial(Mcl1(bf=None, state='M'), Mcl1_0)
    Initial(Noxa(bf=None, state='C'), Noxa_0)
    Initial(CytoC(bf=None, state='M'), CytoC_0)
    Initial(Smac(bf=None, state='M'), Smac_0)

def translocate_tBid_Bax_BclxL():
    """tBid, Bax and BclXL translocate to the mitochondrial membrane."""
    equilibrate(Bid(bf=None, state='T'), Bid(bf=None, state='M'), [1e-1, 1e-3])

    free_Bax = Bax(bf=None, s1=None, s2=None) # Alias for readability
    equilibrate(free_Bax(state='C'), free_Bax(state='M'),
                transloc_rates)

    equilibrate(BclxL(bf=None, state='C'), BclxL(bf=None, state='M'),
                transloc_rates)

def tBid_activates_Bax_and_Bak():
    """tBid activates Bax and Bak."""
    catalyze(Bid(state='M'), Bax(state='M'), Bax(state='A'), bid_effector_rates)
    catalyze(Bid(state='M'), Bak(state='M'), Bak(state='A'), bid_effector_rates)

def tBid_binds_all_anti_apoptotics():
    """tBid binds and inhibits Bcl2, Mcl1, and Bcl-XL."""
    # Doug Green's "MODE 1" inhibition
    bind_table([[                            Bcl2,  BclxL(state='M'),  Mcl1(state='M')],
                [Bid(state='M'),       bcl2_rates,        bcl2_rates,       bcl2_rates]])

def sensitizers_bind_anti_apoptotics():
    """Binding of Bad and Noxa to Bcl2, Mcl1, and Bcl-XL."""
    bind_table([[                       Bcl2,  BclxL(state='M'),  Mcl1(state='M')],
                [Bad(state='M'),  bcl2_rates,        bcl2_rates,             None],
                [Noxa(state='M'),       None,              None,       bcl2_rates]])

def lopez_pore_formation(do_pore_transport=True):
    """ Pore formation and transport process used by all modules.
    """
    alias_model_components()

    # Rates
    pore_max_size = 4
    pore_rates = [[2.040816e-04,  # 1.0e-6/v**2
                   1e-3]] * (pore_max_size - 1)
    pore_transport_rates = [[2.857143e-5, 1e-3, 10]] # 2e-6 / v?

    # Pore formation by effectors
    assemble_pore_sequential(Bax(bf=None, state='A'), pore_max_size, pore_rates)
    assemble_pore_sequential(Bak(bf=None, state='A'), pore_max_size, pore_rates)

    # CytoC, Smac release
    if do_pore_transport:
        pore_transport(Bax(bf=None, state='A'), 4, CytoC(state='M'),
                       CytoC(state='C'), pore_transport_rates)
        pore_transport(Bax(bf=None, state='A'), 4, Smac(state='M'),
                       Smac(state='C'), pore_transport_rates)
        pore_transport(Bak(bf=None, state='A'), 4, CytoC(state='M'),
                       CytoC(state='C'), pore_transport_rates)
        pore_transport(Bak(bf=None, state='A'), 4, Smac(state='M'),
                       Smac(state='C'), pore_transport_rates)

## MOMP Module Implementations ----------------------------------------------

# Embedded Model
def embedded(do_pore_transport=True):
    """ Direct and indirect modes of action, occurring at the membrane.
    """
    alias_model_components()

    declare_initial_conditions()

    translocate_tBid_Bax_BclxL()

    tBid_activates_Bax_and_Bak()

    tBid_binds_all_anti_apoptotics()

    # Autoactivation: Bax and Bak activate their own kind, but only when
    # free (i.e. not part of a pore complex)
    effector_auto_rates = [1e-7, 1e-3, 1]
    catalyze(Bax(active_monomer), Bax(state='M'), Bax(state='A'),
             effector_auto_rates)
    catalyze(Bak(active_monomer), Bak(state='M'), Bak(state='A'),
             effector_auto_rates)

    # Anti-apoptotics bind activated effectors
    # Doug Green's "MODE 2" inhibition
    bind_table([[                            Bcl2,  BclxL(state='M'),        Mcl1],
                [Bax(active_monomer),  bcl2_rates,        bcl2_rates,        None],
                [Bak(active_monomer),        None,        bcl2_rates,  bcl2_rates]])

    sensitizers_bind_anti_apoptotics()

    # Bax and Bak form pores by sequential addition and transport CytoC/Smac
    lopez_pore_formation(do_pore_transport=do_pore_transport)

# Indirect Model
def indirect(do_pore_transport=True):
    """Bax and Bak spontaneously form pores without activation.
       The "activator" tBid binds all of the anti-apoptotics.
    """
    alias_model_components()

    declare_initial_conditions()

    translocate_tBid_Bax_BclxL()

    tBid_binds_all_anti_apoptotics()

    # Bax and Bak spontaneously become activated
    free_Bax = Bax(bf=None, s1=None, s2=None) # Alias
    free_Bak = Bak(bf=None, s1=None, s2=None) # Alias
    equilibrate(free_Bax(state='M'), free_Bax(state='A'), transloc_rates)
    equilibrate(free_Bak(state='M'), free_Bak(state='A'), transloc_rates)

    # Anti-apoptotics bind "active" Bax and Bak to prevent pore formation
    bind_table([[                            Bcl2,  BclxL(state='M'),  Mcl1(state='M')],
                [Bax(active_monomer),  bcl2_rates,        bcl2_rates,             None],
                [Bak(active_monomer),        None,        bcl2_rates,       bcl2_rates]])

    sensitizers_bind_anti_apoptotics()

    # Bax and Bak form pores by sequential addition
    lopez_pore_formation(do_pore_transport=do_pore_transport)

# Direct Model
def direct(do_pore_transport=True):
    """Anti-apoptotics prevent BH3-onlies from activating Bax and Bak.

    Bax and Bak require activation to be able to form pores.
    The anti-apoptotics don't inhibit activated Bax and Bak; their only role
    is to bind BH3-onlies.
    """

    alias_model_components()
    declare_initial_conditions()

    translocate_tBid_Bax_BclxL()

    tBid_activates_Bax_and_Bak()

    tBid_binds_all_anti_apoptotics()

    sensitizers_bind_anti_apoptotics()

    # Bax and Bak form pores by sequential addition
    lopez_pore_formation(do_pore_transport=do_pore_transport)
