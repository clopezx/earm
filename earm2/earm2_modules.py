"""TODO: Write docstring for earm2_modules

More stuff here.
"""

from pysb import *
from pysb.util import alias_model_components
from earm2.macros import *
from pysb.macros import equilibrate
import albeck_modules

# Default site name for binding
site_name = 'bf'

transloc_rates =     [        1e-2, 1e-2]
bcl2_rates =         [1.428571e-05, 1e-3]    # #1.0e-6/v
bid_effector_rates = [        1e-7, 1e-3, 1] # Generalize to catalysis rates?
# Site-value arguments to describe Bax or Bak in the active state but not
# yet oligomerized
active_monomer = {'s1': None, 's2': None, 'state': 'A'}

# MONOMER DECLARATION MACROS ================================================
def all_monomers():
    """Declare the monomers used for the full model, including Bcl-2 proteins.

    Internally calls the macros ligand_to_c8_monomers(),
    bcl2_monomers(), and apaf1_to_parp_monomers() to
    instantiate the monomers for each portion of the pathway.

    The package variable site_name specifies the name of the site to be used
    for all binding reactions (with the exception of Bax and Bak, which have
    additional sites used for oligomerization).

    The 'state' site denotes various localization and/or activity states of a
    Monomer, with 'C' denoting cytoplasmic localization and 'M' mitochondrial
    localization.
    """

    albeck_modules.ligand_to_c8_monomers()
    momp_monomers()
    albeck_modules.apaf1_to_parp_monomers()

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
    # TODO: Bim
    # TODO: Puma
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
    # TODO: Add BclW and Bfl1?
    # == Sensitizers ==================
    Monomer('Bad', [site_name, 'state'], {'state':['C', 'M']})
    Monomer('NOXA', [site_name, 'state'], {'state': ['C', 'M']})
    # TODO: Others???

    # == Cytochrome C and Smac ========
    Monomer('CytoC', [site_name, 'state'], {'state':['M', 'C', 'A']})
    Monomer('Smac', [site_name, 'state'], {'state':['M', 'C', 'A']})

# MOMP SEGMENT ==============================================================

## Macros -------------------------------------------------------------------
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
    pore_transport(Bax, CytoC(state='M'), CytoC(state='C'),
                   pore_transport_rates)
    pore_transport(Bax, Smac(state='M'), Smac(state='C'),
                   pore_transport_rates)
    pore_transport(Bak, CytoC(state='M'), CytoC(state='C'),
                   pore_transport_rates)
    pore_transport(Bak, Smac(state='M'), Smac(state='C'),
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

## MOMP Module Implementations ----------------------------------------------

# Embedded Model
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

# Indirect Model
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

# Direct Model
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
