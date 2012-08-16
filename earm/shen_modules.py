"""PySB implementations of Bcl2-models from the group of Pingping
Shen, along with other derived, closely related models.

In a series of papers from 2007-2010, the research group of Pingping Shen
implemented and investigated models of Bcl-2 family interactions. In this file
we have re-implemented these models using PySB. We have also included a model
from Howells (2011) which is a fairly straightforward extension of a Shen group
model.

Model implementations
---------------------

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
each model.

The models are closely related, and many of the later models are derived from
earlier ones. The models have been implemented in such a way as to make this
hierarchy transparent.

Other functions
---------------

In addition to the implementations of the models themselves, this file also
contains two macros that are re-used by the various models:

- :py:func:`momp_monomers`, which declares the Bcl-2 molecule types used in the
  models, and
- :py:func:`shen_pore_transport`, which declares the set of transport reactions
  required for the release of Cytochrome c and Smac.

"""

from pysb import *
from macros import *
from pysb.macros import catalyze_one_step_reversible, catalyze_one_step, \
                        synthesize_degrade_table, degrade, equilibrate
from pysb.util import alias_model_components

# Default site name for binding
site_name = 'bf'

transloc_rates = [1e-2, 1e-2]

# Useful aliases for Bax state
inactive_monomer = {'state':'C', 's1': None, 's2': None}
active_monomer = {'state':'A', 's1': None, 's2': None}

# For pore_transport sections, drawn from Albeck
v = 0.07
rate_scaling_factor = 1./v

transloc_rates = [1e-2, 1e-2]

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
    # == Sensitizers ==================
    Monomer('Bad',
            [site_name, 'state', 'serine'],
            {'state':['C', 'M'], 'serine':['U', 'P', 'B']})

    # == Cytochrome C and Smac ========
    Monomer('CytoC', [site_name, 'state'], {'state':['M', 'C', 'A']})
    Monomer('Smac', [site_name, 'state'], {'state':['M', 'C', 'A']})

def shen_pore_transport(pore_size=4, micromolar=True):
    """TODO: Document this (hacky) function"""
    if micromolar:
        unit_scaling = 1
    else: # i.e., nanomolar concentrations/rates
        unit_scaling = 1000

    Initial(Smac(state='M', bf=None), Parameter('Smac_0', 0.1*unit_scaling))
    Initial(CytoC(state='M', bf=None), Parameter('CytoC_0', 0.5*unit_scaling))

    # Fwd rate of 2e-6 s^-1 converted to 2 uM^-1 s^-1 (assumes vol of 1.6pL)
    pore_transport(Bax(state='A'), pore_size, Smac(state='M'),
                   Smac(state='C'),
                   [[rate_scaling_factor*2*(1./unit_scaling), 1e-3, 10]])
    pore_transport(Bax(state='A'), pore_size, CytoC(state='M'),
                   CytoC(state='C'),
                   [[rate_scaling_factor*2*(1./unit_scaling), 1e-3, 10]])

## MOMP Module Implementations ---------------------------------------------

def chen_biophys_j(do_pore_assembly=True, do_pore_transport=False):
    """Model drawn from Chen et al. (2007) Biophysical Journal."""
    Parameter('Bcl2_0'  , 1e-1) # Mitochondrial Bcl2
    Parameter('Bax_0'   , 2e-1) # Bax

    alias_model_components()

    Initial(Bax(bf=None, s1=None, s2=None, state='C'), Bax_0)
    Initial(Bcl2(bf=None), Bcl2_0)

    # One-step "kiss-and-run" activation of Bax by tBid
    catalyze_one_step_reversible(
        Bid(state='T', bf=None), Bax(bf=None, **inactive_monomer),
        Bax(bf=None, **active_monomer), [0.5, 1e-1])

    # Bcl2 binds tBid and Bax
    bind_table([[                           Bcl2],
                [Bid(state='T'),       (3, 4e-2)],
                [Bax(active_monomer),  (2, 1e-3)]])

    # Bax can displace Bid from Bcl2
    displace(Bax(active_monomer), Bid(state='T'), Bcl2, 2)

    if do_pore_assembly:
        # Four Bax monomers cooperatively bind to form a tetramer
        assemble_pore_spontaneous(Bax(state='A', bf=None), [2*4, 0])

    if do_pore_transport:
        shen_pore_transport(micromolar=True, pore_size=4)

def chen_febs_indirect(do_pore_assembly=True, do_pore_transport=False):
    """The "indirect" model drawn from Chen et al. (2007) FEBS Letters."""
    Parameter('Bcl2_0'  , 30) # Mitochondrial Bcl2
    Parameter('Bax_0'   , 60) # Bax

    alias_model_components()

    Initial(Bax(bf=None, s1=None, s2=None, state='A'), Bax_0)
    Initial(Bcl2(bf=None), Bcl2_0)

    # (Note: No activation of Bax by tBid, so Bax is in the active state
    # by default)

    # Bcl2 binds tBid and Bax
    bind_table([[                              Bcl2],
                [Bid(state='T'),       (1e-4, 1e-3)],
                [Bax(active_monomer),  (1e-4, 1e-3)]])

    if do_pore_assembly:
        # Four "inactive" Bax monomers cooperatively bind to form a tetramer
        assemble_pore_spontaneous(Bax(state='A', bf=None), [4*1e-3, 1e-3])

    if do_pore_transport:
        shen_pore_transport(micromolar=False, pore_size=4)

def chen_febs_direct(do_pore_assembly=True, do_pore_transport=False):
    """The "direct" model drawn from Chen et al. (2007) FEBS Letters."""
    # All parameters in nanomolar
    # Initial conditions
    Parameter('Bcl2_0' , 30) # Bcl2
    Parameter('Bax_0'  , 60) # InBax
    alias_model_components()
    Initial(Bax(bf=None, s1=None, s2=None, state='C'), Bax_0)
    Initial(Bcl2(bf=None), Bcl2_0)

    # One-step "kiss-and-run" activation of Bax by tBid
    catalyze_one_step_reversible(
        Bid(state='T', bf=None), Bax(bf=None, **inactive_monomer),
        Bax(bf=None, **active_monomer), [1e-3, 1e-3])

    # Bcl2 binds tBid and Bad (a sensitizer) but not Bax
    bind_table([[                         Bcl2],
                [Bid(state='T'),  (1e-4, 1e-3)],
                [Bad(state='M'),  (1e-4, 1e-3)]])

    if do_pore_assembly:
        # Four Bax monomers cooperatively bind to form a tetramer
        assemble_pore_spontaneous(Bax(state='A', bf=None), [4*1e-3, 1e-3])

    if do_pore_transport:
        shen_pore_transport(micromolar=False, pore_size=4)

def cui_direct(do_pore_transport=False):
    """The "direct" model drawn from Cui et al. (2008) PLoS One."""
    # All parameters in nanomolar
    # Build on the direct model from Chen et al. (2007) FEBS Lett. by:
    chen_febs_direct(do_pore_assembly=False,
                        do_pore_transport=False)
    alias_model_components()

    # 1. Overriding some parameter values,
    one_step_BidT_BaxC_to_BidT_BaxA_kf.value = 0.0005 # originally 0.001
    bind_BidT_Bcl2_kf.value = 0.001                   # originally 0.0001

    # 2. Adding a Bad-for-Bid displacement reaction,
    displace_reversibly(Bad(state='M'), Bid(state='T'), Bcl2,
                        [0.0001, 0.001])

    # 3. Adding simplified MAC formation (Bax dimerization)
    # NOTE: Even though this binding reaction is homomeric (which would
    # imply that the forward rate should be divided by two) the fact that the
    # binding reaction can occur in two different ways (s1 on one Bax binding
    # to the s2 on another, vs. s2 on the first Bax binding to s1 on the second)
    # requires that the rate be scaled back by a factor of two. These two
    # scaling factors cancel out, so the forward rate constant is not scaled,
    # and is used in the ODE with its nominal value.
    active_unbound = {'state': 'A', 'bf': None}
    assemble_pore_sequential(Bax(**active_unbound), 2, [[0.0002, 0.02]])

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
        shen_pore_transport(pore_size=2, micromolar=False)

def cui_direct1(do_pore_transport=False):
    """The "direct 1" model drawn from Cui et al. (2008) PLoS One."""
    alias_model_components()

    # Build on the base "direct" model...
    cui_direct(do_pore_transport=do_pore_transport)

    # ...by adding inhibition of Bax by Bcl2,
    bind(Bax(state='A', s1=None, s2=None), Bcl2, [0.005, 0.001])

    # ...associated displacement reactions
    displace_reversibly(Bax(active_monomer), Bid(state='T'), Bcl2,
                        [0.005, 0.001])
    displace_reversibly(Bad(state='M'), Bax(active_monomer), Bcl2,
                        [0.0001, 0.005])

    # ...and degradation of the active Bax:Bcl2 complex
    degrade(Bax(bf=1) % Bcl2(bf=1), 0.005)

def cui_direct2(do_pore_transport=False):
    """The "direct 2" model drawn from Cui et al. (2008) PLoS One."""
    alias_model_components()

    # Build on the "direct 1" model...
    cui_direct1(do_pore_transport=do_pore_transport)

    # By adding simultaneous auto-activation and dimerization of Bax
    Rule('Bax_autoactivation_dimerization',
        Bax(state='A', bf=None, s1=None, s2=None) +
        Bax(state='C', bf=None, s1=None, s2=None) >>
        Bax(state='A', bf=None, s1=1, s2=None) %
        Bax(state='A', bf=None, s1=None, s2=1),
        Parameter('Bax_autoactivation_dimerization_k', 0.0002))

def howells(do_pore_assembly=True, do_pore_transport=False):
    """The model drawn from Howells et al. (2011) J. Theor. Biol."""
    # Build on the model from Chen et al. (2007) Biophys J:
    chen_biophys_j(do_pore_assembly=do_pore_assembly,
                     do_pore_transport=do_pore_transport)
    alias_model_components()

    # Override a few parameter values from the pre-existing model
    bind_BidT_Bcl2_kr.value = 2e-3 # was 4e-2 in Chen 2007 Biophys J
    bind_BaxA_Bcl2_kr.value = 2e-3 # was 1e-3 in Chen 2007 Biophys J
    spontaneous_pore_BaxA_to_Bax4_kf.value = 2000*4 # was 2 in Chen 2007 B.J.
    spontaneous_pore_BaxA_to_Bax4_kr.value = 5e-5   # was 0 in Chen 2007 B.J.

    # Translocation equilibrium between unphosphorylated cytosolic and
    # mitochondrial Bad
    equilibrate(Bad(state='C', serine='U', bf=None),
                Bad(state='M', serine='U', bf=None), [1e-2, 2e-3])

    # Bad binds Bcl2
    bind(Bad(state='M'), Bcl2, [15, 2e-3])

    # Bad displaces tBid from Bcl2
    displace(Bad(state='M'), Bid(state='T'), Bcl2, 5) # k_tBid_rel1 in paper

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

    # Sequester phospho-Bad by "binding" 14-3-3 domains
    Rule('sequester_BadCP_to_BadC1433',
         Bad(state='C', serine='P') >> Bad(state='C', serine='B'),
         Parameter('sequester_BadCP_to_BadC1433_k', 1e-3))

    # Release of Bad from 14-3-3 domains
    Rule('release_BadC1433_to_BadCU',
         Bad(state='C', serine='B') >> Bad(state='C', serine='U'),
         Parameter('release_BadC1433_to_BadCU_k', 8.7e-4))
