"""
Overview
========

Three models of MOMP (:py:func:`direct`, :py:func:`indirect`, and
:py:func:`embedded`), each incorporating a larger repertoire of Bcl-2 family
members than previously published models, including:

* One **activator,** Bid.
* Two **sensitizers,** Bad and Noxa.
* Two **effectors,** Bax and Bak.
* Three **anti-apoptotics**, Bcl-2, Bcl-xL, and Mcl-1.

The Models
----------

Note that in each of the three models, interactions between Bcl-2 proteins only
occur at the mitochondrial membrane. The following are brief descriptions of
each model.

* :py:func:`direct`. In this model, tBid directly activates both Bax and Bak;
  the anti-apoptotics bind tBid and the sensitizers (Bad and Noxa) but not
  Bax and Bak.
* :py:func:`indirect`. Bax and Bak are not explicitly activated by tBid, but
  rather are in an equilibrium between inactive and active states. The
  anti-apoptotics bind tBid, sensitizers, and Bax and Bak.
* :py:func:`embedded`. Combines elements of both direct and indirect: tBid
  activates Bax and Bak; the anti-apoptotics bind tBid, sensitizers and Bax and
  Bak. In addition, Bax and Bak are able to auto-activate.

Organization of models into Motifs
----------------------------------

Because the three models share many aspects, the mechanisms that they share have
been written as small "motifs" implemented as subroutines. These are:

* :py:func:`translocate_tBid_Bax_BclxL`
* :py:func:`tBid_activates_Bax_and_Bak`
* :py:func:`tBid_binds_all_anti_apoptotics`
* :py:func:`sensitizers_bind_anti_apoptotics`
* :py:func:`effectors_bind_anti_apoptotics`
* :py:func:`lopez_pore_formation`

The implementation details of these motifs can be seen by examining the
source code.

Monomer and initial declaration functions
-----------------------------------------

The models share the same set of Monomer and initial condition declarations,
which are contained with the following two functions:

* :py:func:`momp_monomers`
* :py:func:`declare_initial_conditions`
"""

# Preliminaries
# =============
#
# We'll need everything from the pysb core, and some macros:

from pysb import *
from shared import *
from pysb.macros import equilibrate
from pysb.util import alias_model_components

# Globals
# -------

# A variety of default rate constant values
bcl2_rates =       [1.428571e-05, 1e-3]    # 1.0e-6/v
activation_rates = [        1e-7, 1e-3, 1] 

# Shared functions
# ================

# Monomer and initial condition declarations
# ------------------------------------------

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

    # **Activators**.
    # Bid, states: Untruncated, Truncated, truncated and Mitochondrial
    Monomer('Bid', [site_name, 'state'], {'state':['U', 'T', 'M']})
    # **Effectors**
    # Bax, states: Cytoplasmic, Mitochondrial, Active
    # sites 's1' and 's2' are used for pore formation
    Monomer('Bax', [site_name, 's1', 's2', 'state'], {'state':['C', 'M', 'A']})
    # Bak, states: inactive and Mitochondrial, Active (and mitochondrial)
    # sites 's1' and 's2' are used for pore formation
    Monomer('Bak', [site_name, 's1', 's2', 'state'], {'state':['M', 'A']})
    # **Anti-Apoptotics**
    Monomer('Bcl2', [site_name])
    Monomer('BclxL', [site_name, 'state'], {'state':['C', 'M']})
    Monomer('Mcl1', [site_name, 'state'], {'state':['M', 'C']})
    # **Sensitizers**
    Monomer('Bad', [site_name, 'state'], {'state':['C', 'M']})
    Monomer('Noxa', [site_name, 'state'], {'state': ['C', 'M']})

    # **Cytochrome C and Smac**
    Monomer('CytoC', [site_name, 'state'], {'state':['M', 'C', 'A']})
    Monomer('Smac', [site_name, 'state'], {'state':['M', 'C', 'A']})

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

# Motifs
# ------

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
    catalyze(Bid(state='M'), Bax(state='M'), Bax(state='A'), activation_rates)
    catalyze(Bid(state='M'), Bak(state='M'), Bak(state='A'), activation_rates)

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

def effectors_bind_anti_apoptotics():
    """Binding of Bax and Bak to Bcl2, BclxL, and Mcl1."""
    bind_table([[                            Bcl2,  BclxL(state='M'),        Mcl1],
                [Bax(active_monomer),  bcl2_rates,        bcl2_rates,        None],
                [Bak(active_monomer),        None,        bcl2_rates,  bcl2_rates]])

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

# MOMP model implementations
# ==========================

def embedded(do_pore_transport=True):
    """ Direct and indirect modes of action, occurring at the membrane.
    """
    alias_model_components()

    declare_initial_conditions()

    translocate_tBid_Bax_BclxL()

    tBid_activates_Bax_and_Bak()


    # Autoactivation: Bax and Bak activate their own kind, but only when
    # free (i.e. not part of a pore complex)
    catalyze(Bax(active_monomer), Bax(state='M'), Bax(state='A'),
             activation_rates)
    catalyze(Bak(active_monomer), Bak(state='M'), Bak(state='A'),
             activation_rates)

    # Anti-apoptotics bind activator tBid
    # Doug Green's "MODE 1" inhibition
    tBid_binds_all_anti_apoptotics()

    # Anti-apoptotics bind activated effectors
    # Doug Green's "MODE 2" inhibition
    effectors_bind_anti_apoptotics()

    sensitizers_bind_anti_apoptotics()

    # Bax and Bak form pores by sequential addition and transport CytoC/Smac
    lopez_pore_formation(do_pore_transport=do_pore_transport)

def indirect(do_pore_transport=True):
    """Bax and Bak spontaneously form pores without activation.
       The "activator" tBid binds all of the anti-apoptotics.
    """
    alias_model_components()

    declare_initial_conditions()

    translocate_tBid_Bax_BclxL()

    # Bax and Bak spontaneously become activated
    free_Bax = Bax(bf=None, s1=None, s2=None) # Alias
    free_Bak = Bak(bf=None, s1=None, s2=None) # Alias
    equilibrate(free_Bax(state='M'), free_Bax(state='A'), transloc_rates)
    equilibrate(free_Bak(state='M'), free_Bak(state='A'), transloc_rates)

    # Anti-apoptotics bind activator tBid
    # Doug Green's "MODE 1" inhibition
    tBid_binds_all_anti_apoptotics()

    # Anti-apoptotics bind activated effectors
    # Doug Green's "MODE 2" inhibition
    effectors_bind_anti_apoptotics()

    sensitizers_bind_anti_apoptotics()

    # Bax and Bak form pores by sequential addition
    lopez_pore_formation(do_pore_transport=do_pore_transport)

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

    # Anti-apoptotics bind activator tBid
    # Doug Green's "MODE 1" inhibition
    tBid_binds_all_anti_apoptotics()

    sensitizers_bind_anti_apoptotics()

    # Bax and Bak form pores by sequential addition
    lopez_pore_formation(do_pore_transport=do_pore_transport)
