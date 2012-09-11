"""
Overview
========

PySB implementations of Bcl2-models from the group of Pingping Shen, along with
other derived, closely related models.

In a series of papers from 2007-2010, the research group of Pingping Shen
implemented and investigated models of Bcl-2 family interactions. In this file
we have re-implemented these models using PySB. We have also included a model
from [Howells2011]_ which is a fairly straightforward extension of a
Shen group model from [Chen2007biophysj]_.

MOMP model implementations
--------------------------

The implementations of the various models are contained within the
following functions:

- :py:func:`chen_biophys_j`
- :py:func:`chen_febs_indirect`
- :py:func:`chen_febs_direct`
- :py:func:`cui_direct`
- :py:func:`cui_direct1`
- :py:func:`cui_direct2`
- :py:func:`howells`

Model descriptions (with references) are available in the documentation for
each function.

The models are closely related, and many of the later models are derived from
earlier ones. The models have been implemented in such a way as to make this
hierarchy transparent.

Shared functions
----------------

In addition to the implementations of the models themselves, this file also
contains two macros that are re-used by the various models:

- :py:func:`momp_monomers`, which declares the Bcl-2 molecule types used in the
  models, and
- :py:func:`shen_pore_transport`, which declares the set of transport reactions
  required for the release of Cytochrome c and Smac.

Parameter values
----------------

In the original papers, species quantities and forward rate constants were
either given in units of micromolar ([Chen2007biophysj]_, [Howells2011]_) or
nanomolar ([Chen2007febs]_, [Cui2008]_).  For consistency, these values have
been written in terms of their molar equivalents (for example, 0.1uM = 0.1e-6 M;
3 uM^-1 s^1 = 3e6 M^-1 s^-1).  Concentrations have been converted into units of
numbers of molecules according to:

No. of molecules = Conc * N_A * vol

where N_A is Avogadro's number and vol is the cell volume, which is given a
default value in the global variable `V` defined in :doc:`shared`.
Similarly, forward rate constants are converted into stochastic rate constants
according to:

Stoch. rate constant = Det. rate constant / (N_A * vol)
"""

# Preliminaries
# =============

# We'll need everything from the pysb core and some macros:

from pysb import *
from shared import *
from pysb.macros import catalyze_one_step_reversible, catalyze_one_step, \
                        synthesize_degrade_table, degrade, equilibrate
from pysb.util import alias_model_components

# Shared functions
# ================

# The Shen model functions share some logic that is contained in the following
# two functions.

def momp_monomers():
    """Declares the signatures of the Bcl-2 family monomers used in all of the
    Shen models.

    In principle, each Shen MOMP model implementation could declare its own
    set of Bcl-2 monomers, each with its own site and state signature. In the
    interest of consistency, a unified set of monomer signatures that supports
    all of the models is defined here.
    """

    # **Activators.** We use Bid as a representative for the generic
    # "activator" "Act" or "BH3" used in the Shen models. Bid has three states:
    # Untruncated, Truncated, and truncated and Mitochondrial.
    Monomer('Bid', ['bf', 'state'], {'state':['U', 'T', 'M']})

    # **Effector***. Bax, states: Cytoplasmic, Mitochondrial, Active.
    # Sites 's1' and 's2' are used for pore formation.
    Monomer('Bax', ['bf', 's1', 's2', 'state'], {'state':['C', 'M', 'A']})

    # **Anti-apoptotic**. Bcl-2 is considered to be constitutively
    # mitochondrial so it doesn't have a 'state' site, only the single binding
    # site.
    Monomer('Bcl2', ['bf'])

    # **Sensitizer**. We use Bad as a representative for the generic "Ena"
    # ("Enabler," another term for sensitizer) used in the Shen models. Bad can
    # be both Cytosolic and Mitochondrial. To support the Bad phosphorylation
    # model of [Howells2011]_ we also add a site "serine", with three
    # states (Unphosphorylated, Phosphorylated, and Bound (i.e., to 14-3-3
    # domains).  These are only used in the Howells model--for all other
    # models, the serine is set as unphosphorylated (i.e., serine='U').
    Monomer('Bad',
            ['bf', 'state', 'serine'],
            {'state':['C', 'M'], 'serine':['U', 'P', 'B']})

    # **Cytochrome  C and Smac**. We include these species here so that the
    # the Shen models can optionally implement CytoC/Smac release. This allows
    # them to interface with the downstream apoptotic machinery when composed
    # into a larger model. Both Cytochrome C and Smac have three states:
    # Mitochondrial (the initial state), Cytosolic (for after release) and
    # Active.
    Monomer('CytoC', ['bf', 'state'], {'state':['M', 'C', 'A']})
    Monomer('Smac', ['bf', 'state'], {'state':['M', 'C', 'A']})

def shen_pore_transport(pore_size=4):
    # TODO: Move this to albeck_modules?
    """Implements release of Cytochrome C and Smac.

    Uses the same model as the original EARM 1.0 ([Albeck2008]_), in
    which pore transport is modeled as binding of the cargo
    (cytochrome C or Smac) to the active pore, and then release, in a
    catalysis-like mechanism.

    The initial conditions for cytochrome C and Smac, and the rate constants
    for transport, are also taken from EARM 1.0.
    """

    Initial(Smac(state='M', bf=None), Parameter('Smac_0', 1e5))
    Initial(CytoC(state='M', bf=None), Parameter('CytoC_0', 5e5))

    pore_transport(Bax(state='A'), pore_size, Smac(state='M'),
                   Smac(state='C'),
                   [[rate_scaling_factor*2e-6, 1e-3, 10]])
    pore_transport(Bax(state='A'), pore_size, CytoC(state='M'),
                   CytoC(state='C'),
                   [[rate_scaling_factor*2e-6, 1e-3, 10]])

# MOMP model implementations
# ==========================

def chen_biophys_j(do_pore_assembly=True, do_pore_transport=False):
    """Model drawn from [Chen2007biophysj]_.

    Model features (see the source code):

    * Activation of Bax by an activator (tBid) in a one-step, hit-and-run
      manner; Bax activation is reversible.
    * Bcl2 binds both tBid and Bax Bax can displace tBid from Bcl-2 (but not
      the reverse).
    * If Bax oligomerization is incorporated into the model (see
      `do_pore_assembly` argument, below), then this occurs as a spontaneous,
      order 4 reaction.

    This model combines both "direct" type and "indirect" type elements in
    that Bcl-2 is capable of binding both Bid and Bax (see `bind_table` call
    in the source code).

    Parameters
    ----------
    do_pore_assembly : True (default) or False
        If True, adds the formation of Bax oligomers to the model. If False,
        the model's most downstream element is Bax activation. This is included
        for two reasons: first, the original publication included two variant
        models, one with and one without Bax oligomerization, so this allows
        this aspect of the original models to be explored. Second, it allows
        a model that extends this model to implement a different model of Bax
        pore assembly (for example, as is the case with cui_direct).
    do_pore_transport : True or False (default)
        If True, adds the release of Cytochrome C and Smac to the model by
        calling the function :py:func:`shen_pore_transport`. If CytoC/Smac
        release are not incorporated into the model, the model matches the
        originally published model but can't be composed into the full
        extrinsic apoptosis pathway.
    """

    Parameter('Bcl2_0', 0.1e-6 * N_A * V) # Mitochondrial Bcl2
    Parameter('Bax_0',  0.2e-6 * N_A * V) # Bax

    alias_model_components()

    # Bax is in the Cytosolic, inactive state by default.
    Initial(Bax(bf=None, s1=None, s2=None, state='C'), Bax_0)
    Initial(Bcl2(bf=None), Bcl2_0)

    # One-step "kiss-and-run" activation of Bax by tBid:
    catalyze_one_step_reversible(
        Bid(state='T', bf=None), Bax(bf=None, **inactive_monomer),
        Bax(bf=None, **active_monomer), [0.5e6/(N_A*V), 1e-1])

    # Bcl2 binds tBid and Bax:
    bind_table([[                                            Bcl2],
                [Bid(state='T'),       (3e6/(N_A*V), 4e-2)],
                [Bax(active_monomer),  (2e6/(N_A*V), 1e-3)]])

    # Bax can displace Bid from Bcl2:
    displace(Bax(active_monomer), Bid(state='T'), Bcl2, 2e6/(N_A*V))

    if do_pore_assembly:
        # Four Bax monomers cooperatively bind to form a tetramer
        assemble_pore_spontaneous(Bax(state='A', bf=None),
                                  [2e6*4/(N_A*V), 0])

    if do_pore_transport:
        # Release Cytochrome C and Smac:
        shen_pore_transport(pore_size=4)

def chen_febs_indirect(do_pore_assembly=True, do_pore_transport=False):
    """The "indirect" model drawn from [Chen2007febs]_.

    Model features (see the source code):

    * There is no activation of Bax by tBid. Bax starts out constitutively
      "active" in that in its initial state, it is able to form oligomers.
    * Bcl-2 can bind tBid and Bax.

    Parameters
    ----------
    do_pore_assembly : True (default) or False
        As for :py:func:`chen_biophys_j`.
    do_pore_transport : True or False (default)
        As for :py:func:`chen_biophys_j`.
    """

    Parameter('Bcl2_0'  , 30e-9) # Mitochondrial Bcl2
    Parameter('Bax_0'   , 60e-9) # Bax

    alias_model_components()

    Initial(Bax(bf=None, s1=None, s2=None, state='A'), Bax_0)
    Initial(Bcl2(bf=None), Bcl2_0)

    # (Note: No activation of Bax by tBid, so Bax is in the active state
    # by default)

    # Bcl2 binds tBid and Bax:
    bind_table([[                                            Bcl2],
                [Bid(state='T'),       (1e5/(N_A*V), 1e-3)],
                [Bax(active_monomer),  (1e5/(N_A*V), 1e-3)]])

    if do_pore_assembly:
        # Four "inactive" Bax monomers cooperatively bind to form a tetramer:
        assemble_pore_spontaneous(Bax(state='A', bf=None),
                                  [4*1e6/(N_A*V), 1e-3])

    if do_pore_transport:
        # Release Cytochrome C and Smac:
        shen_pore_transport(pore_size=4)

def chen_febs_direct(do_pore_assembly=True, do_pore_transport=False):
    """The "direct" model drawn from [Chen2007febs]_.

    Model features (see the source code):

    * Activation of Bax by an activator (tBid) in a one-step, hit-and-run
      manner; Bax activation is reversible.
    * Bcl-2 can bind tBid, but not Bax.

    Parameters
    ----------
    do_pore_assembly : True (default) or False
        As for :py:func:`chen_biophys_j`.
    do_pore_transport : True or False (default)
        As for :py:func:`chen_biophys_j`.
    """

    # Initial conditions
    Parameter('Bcl2_0' , 30e-9) # Bcl2
    Parameter('Bax_0'  , 60e-9) # InBax

    alias_model_components()

    Initial(Bax(bf=None, s1=None, s2=None, state='C'), Bax_0)
    Initial(Bcl2(bf=None), Bcl2_0)

    # One-step "kiss-and-run" activation of Bax by tBid
    catalyze_one_step_reversible(
        Bid(state='T', bf=None), Bax(bf=None, **inactive_monomer),
        Bax(bf=None, **active_monomer), [1e6/(N_A*V), 1e-3])

    # Bcl2 binds tBid and Bad (a sensitizer) but not Bax
    bind_table([[                                       Bcl2],
                [Bid(state='T'),  (1e5/(N_A*V), 1e-3)],
                [Bad(state='M'),  (1e5/(N_A*V), 1e-3)]])

    if do_pore_assembly:
        # Four Bax monomers cooperatively bind to form a tetramer
        assemble_pore_spontaneous(Bax(state='A', bf=None),
                                  [4*1e6/(N_A*V), 1e-3])

    if do_pore_transport:
        # Release Cytochrome C and Smac:
        shen_pore_transport(pore_size=4)

def cui_direct(do_pore_transport=False):
    """The "direct" model drawn from [Cui2008]_.

    Builds on the direct model from [Chen2007febs]_, implemented
    in :py:func:`chen_febs_direct` (see source code).

    Parameters
    ----------
    do_pore_transport : True or False (default)
        As for :py:func:`chen_biophys_j`.
    """

    # Build on the direct model from [Chen2007febs]_ by...
    # (NOTE that we set the keyword argument `do_pore_assembly` to False here:
    # this is because we don't want to use a tetrameric pore as was implemented
    # in the Chen model, we want to implement a dimeric pore. See comment #3,
    # below.)
    chen_febs_direct(do_pore_assembly=False, do_pore_transport=False)

    alias_model_components()

    # 1. Overriding some parameter values:
    # Cut the activation rate of Bax by tBid by half (was originally 1e6):
    one_step_BidT_BaxC_to_BidT_BaxA_kf.value = 5e5/(N_A*V)
    # Make Bid/Bcl2 binding 10-fold tighter (was originally 1e5)
    bind_BidT_Bcl2_kf.value = 1e6/(N_A*V)

    # 2. Adding a Bad-for-Bid displacement reaction,
    displace_reversibly(Bad(state='M'), Bid(state='T'), Bcl2,
                        [1e5/(N_A*V), 0.001])

    # 3. Adding simplified MAC formation (Bax dimerization)
    # NOTE: Even though this binding reaction is homomeric (which would
    # imply that the forward rate should be divided by two) the fact that the
    # binding reaction can occur in two different ways (s1 on one Bax binding
    # to the s2 on another, vs. s2 on the first Bax binding to s1 on the second)
    # requires that the rate be scaled back by a factor of two. These two
    # scaling factors cancel out, so the forward rate constant is not scaled,
    # and is used in the ODE with its nominal value.
    active_unbound = {'state': 'A', 'bf': None}
    assemble_pore_sequential(Bax(**active_unbound), 2,
                                 [[2e5/(N_A*V), 0.02]])

    # 4. Adding synthesis and degradation reactions
    Bax2 = Bax(s1=1, s2=None) % Bax(s1=None, s2=1)
    synthesize_degrade_table(
       [[Bax(bf=None, **inactive_monomer),      0.06,  0.001],
        [Bax(bf=None, **active_monomer),        None,  0.001],
        [Bid(state='T', bf=None),              0.001,  0.001],
        [Bcl2(bf=None),                         0.03,  0.001],
        [Bid(state='T', bf=1) % Bcl2(bf=1),     None,  0.005],
        [Bad(state='M', bf=None, serine='U'),  0.001,  0.001],
        [Bad(bf=1) % Bcl2(bf=1),                None,  0.005],
        [Bax2,                                  None, 0.0005]])

    if do_pore_transport:
        # Release Cytochrome C and Smac:
        shen_pore_transport(pore_size=2)

def cui_direct1(do_pore_transport=False):
    """The "direct 1" model drawn from [Cui2008]_.

    Builds on the (base) direct model from [Cui2008]_,
    implemented in :py:func:`cui_direct` (see source code).

    Parameters
    ----------
    do_pore_transport : True or False (default)
        As for :py:func:`chen_biophys_j`.
    """

    alias_model_components()

    # Build on the base "direct" model...
    cui_direct(do_pore_transport=do_pore_transport)

    # ...by adding inhibition of Bax by Bcl2,
    bind(Bax(state='A', s1=None, s2=None), Bcl2, [5e6/(N_A*V), 0.001])

    # ...associated displacement reactions
    displace_reversibly(Bax(active_monomer), Bid(state='T'), Bcl2,
                        [5e6/(N_A*V), 0.001])
    displace_reversibly(Bad(state='M'), Bax(active_monomer), Bcl2,
                        [1e5/(N_A*V), 0.005])

    # ...and degradation of the active Bax:Bcl2 complex
    degrade(Bax(bf=1) % Bcl2(bf=1), 0.005)

def cui_direct2(do_pore_transport=False):
    """The "direct 2" model drawn from [Cui2008]_.

    Builds on the "direct 1" model from [Cui2008]_, implemented
    in :py:func:`cui_direct1` (see source code).

    Parameters
    ----------
    do_pore_transport : True or False (default)
        As for :py:func:`chen_biophys_j`.
    """
    alias_model_components()

    # Build on the "direct 1" model...
    cui_direct1(do_pore_transport=do_pore_transport)

    # By adding simultaneous auto-activation and dimerization of Bax
    Rule('Bax_autoactivation_dimerization',
        Bax(state='A', bf=None, s1=None, s2=None) +
        Bax(state='C', bf=None, s1=None, s2=None) >>
        Bax(state='A', bf=None, s1=1, s2=None) %
        Bax(state='A', bf=None, s1=None, s2=1),
        Parameter('Bax_autoactivation_dimerization_k', 2e5/(N_A*V)))

def howells(do_pore_assembly=True, do_pore_transport=False):
    """The model drawn from [Howells2011]_.

    This model builds on the model from [Chen2007biophysj]_,
    implemented in :py:func:`chen_biophys_j`. The core reactions from the Chen
    et al. model are the same, but Howells et al. modify some parameter values
    and add a number of Bad-related reactions, including (see source code):

    * Unphosphorylated Bad spontaneously translocates between cytosol and
      mitochondria
    * Bad binds Bcl-2
    * Bad displaces tBid from Bcl-2
    * Cytosolic, mitochondrial, and Bad in a mitochondrial Bad:Bcl2 complex
      can be phosphorylated at various rates (this is modeled as a first-order
      reaction with no explicit representation of kinases)
    * Bad can be sequestered by, and released from, 14-3-3 domains in the
      cytosol (modeled as a first-order reaction with no explicit
      representation of 14-3-3-domain-containing proteins)

    Parameters
    ----------
    do_pore_assembly : True (default) or False
        As for :py:func:`chen_biophys_j`.
    do_pore_transport : True or False (default)
        As for :py:func:`chen_biophys_j`.
    """

    # Build on the model from [Chen2007biophysj]_:
    chen_biophys_j(do_pore_assembly=do_pore_assembly,
                     do_pore_transport=do_pore_transport)
    alias_model_components()

    # Override a few parameter values from the pre-existing model
    bind_BidT_Bcl2_kr.value = 2e-3 # was 4e-2 in Chen 2007 Biophys J
    bind_BaxA_Bcl2_kr.value = 2e-3 # was 1e-3 in Chen 2007 Biophys J
    spontaneous_pore_BaxA_to_Bax4_kf.value = 2000e6*4/(N_A*V)
                                                    # was 2e6 in Chen 2007 B.J.
    spontaneous_pore_BaxA_to_Bax4_kr.value = 5e-5   # was 0 in Chen 2007 B.J.

    # Translocation equilibrium between unphosphorylated cytosolic and
    # mitochondrial Bad
    equilibrate(Bad(state='C', serine='U', bf=None),
                Bad(state='M', serine='U', bf=None), [1e-2, 2e-3])

    # Bad binds Bcl2
    bind(Bad(state='M'), Bcl2, [15e6/(N_A*V), 2e-3])

    # Bad displaces tBid from Bcl2 (parameter `k_tBid_rel1` in paper)
    displace(Bad(state='M'), Bid(state='T'), Bcl2, 5e6/(N_A*V))

    # Phosphorylation of Bad
    phosphorylate_Bad_k1 = Parameter('phosphorylate_Bad_k1', 1e-3)
    phosphorylate_Bad_k2 = Parameter('phosphorylate_Bad_k2', 1e-4)
    Rule('phosphorylate_BadCU_to_BadCP',     # Cytosolic Bad
         Bad(state='C', serine='U') >> Bad(state='C', serine='P'),
         phosphorylate_Bad_k1)
    Rule('phosphorylate_BadMU_to_BadCP',     # Mitochondrial Bad
         Bad(state='M', serine='U', bf=None) >>
         Bad(state='C', serine='P', bf=None),
         phosphorylate_Bad_k1)
    Rule('phosphorylate_BadMUBcl2_to_BadCP', # Mitochondrial Bad:Bcl2
         Bad(state='M', serine='U', bf=1) % Bcl2(bf=1) >>
         Bad(state='C', serine='P', bf=None) + Bcl2(bf=None),
         phosphorylate_Bad_k2)

    # Sequester phospho-Bad by "binding" 14-3-3 domains (parameter `k_BAD_seq`)
    Rule('sequester_BadCP_to_BadC1433',
         Bad(state='C', serine='P') >> Bad(state='C', serine='B'),
         Parameter('sequester_BadCP_to_BadC1433_k', 1e-3))

    # Release of Bad from 14-3-3 domains (parameter `k_BAD_rel`)
    Rule('release_BadC1433_to_BadCU',
         Bad(state='C', serine='B') >> Bad(state='C', serine='U'),
         Parameter('release_BadC1433_to_BadCU_k', 8.7e-4))


