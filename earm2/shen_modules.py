from pysb import *
from macros import *
from pysb.macros import catalyze_one_step_reversible, synthesize_degrade_table
from pysb.util import alias_model_components
#from pysb.macros import catalyze_one_step_reversible

# Nomenclature from the paper:
# Ena = Bad
# Act = tBid
# Bcl2 = Bcl2
# InBax = Bax(state='C')
# AcBax = Bax(state='A')
# MAC = Pore()

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
    Monomer('Bad', [site_name, 'state'], {'state':['C', 'M']})

    # == Cytochrome C and Smac ========
    Monomer('CytoC', [site_name, 'state'], {'state':['M', 'C', 'A']})
    Monomer('Smac', [site_name, 'state'], {'state':['M', 'C', 'A']})

def momp_initial_conditions(bid_state='U'):
    # TODO: Eliminate this function?
    """Declare initial conditions for Bcl-2 family proteins, Cyto c, and Smac.

    Parameters
    ----------
    bid_state : string; 'U' or 'T'
        Specifies the initial state of Bid. 'U' indicates un-truncated (i.e.,
        full-length); 'T' indicates truncated.
    """

    Parameter('Bid_0'   , 4.0e4) # Bid
    Parameter('Bcl2_0'  , 2.0e4) # Mitochondrial Bcl2
    Parameter('Bad_0'   , 1.0e3) # Bad
    Parameter('Bax_0'   , 0.8e5) # Bax
    Parameter('CytoC_0' , 5.0e5) # cytochrome c
    Parameter('Smac_0'  , 1.0e5) # Smac

    alias_model_components()

    Initial(Bid(bf=None, state=bid_state), Bid_0)
    Initial(Bad(bf=None, state='C'), Bad_0)
    Initial(Bax(bf=None, s1=None, s2=None, state='C'), Bax_0)
    Initial(Bcl2(bf=None), Bcl2_0)
    Initial(CytoC(bf=None, state='M'), CytoC_0)
    Initial(Smac(bf=None, state='M'), Smac_0)

def chen2007BiophysJ(pore_assembly=False, pore_transport=False):
    # TODO: change all initial conditions and param values to Molar
    Parameter('Bid_0'   , 0) # Bid
    Parameter('Bcl2_0'  , 1e-1) # Mitochondrial Bcl2
    Parameter('Bax_0'   , 2e-1) # Bax

    alias_model_components()

    Initial(Bid(bf=None, state='T'), Bid_0)
    Initial(Bax(bf=None, s1=None, s2=None, state='C'), Bax_0)
    Initial(Bcl2(bf=None), Bcl2_0)

    # Aliases for Bax state
    inactive_monomer = {'state':'C', 's1': None, 's2': None}
    active_monomer = {'state':'A', 's1': None, 's2': None}

    # One-step "kiss-and-run" activation of Bax by tBid
    catalyze_one_step_reversible(
        Bid(state='T', bf=None), Bax(bf=None, **inactive_monomer),
        Bax(bf=None, **active_monomer), [0.5, 1e-1])

    # Bcl2 binds tBid and Bax
    bind_table([[                              Bcl2],
                [Bid(state='T'),       (3, 4e-2)],
                [Bax(active_monomer),  (2, 1e-3)]]) # was 1e-3 for kr

    # Bax can displace Bid from Bcl2
    displace_reversibly(Bax(active_monomer), Bid(state='T'), Bcl2, [2, 0])

    if pore_assembly:
        # Four Bax monomers cooperatively bind to form a tetramer
        assemble_pore_spontaneous(Bax(state='A', bf=None), [2, 0])

def chen2007FEBS_indirect(pore_assembly=False):
    # TODO: change all initial conditions and param values to Molar
    Parameter('Bid_0'   , 0) # Bid
    Parameter('Bcl2_0'  , 30) # Mitochondrial Bcl2
    Parameter('Bax_0'   , 60) # Bax

    alias_model_components()

    Initial(Bid(bf=None, state='T'), Bid_0)
    Initial(Bax(bf=None, s1=None, s2=None, state='C'), Bax_0)
    Initial(Bcl2(bf=None), Bcl2_0)

    # (Note: No activation of Bax by tBid, so binding and oligomerization
    # reactions apply to the unactivated, or state='C' form of Bax)

    # Bcl2 binds tBid and Bax
    monomeric = {'state': 'C', 's1': None, 's2': None}
    bind_table([[                         Bcl2],
                [Bid(state='T'),  (1e-4, 1e-3)],
                [Bax(monomeric),  (1e-4, 1e-3)]])

    if pore_assembly:
        # Four "inactive" Bax monomers cooperatively bind to form a tetramer
        assemble_pore_spontaneous(Bax(state='C', bf=None), [4*1e-3, 1e-3])

def chen2007FEBS_direct(pore_assembly=False):
    # TODO: change all initial conditions and param values to Molar
    Parameter('Bid_0'   , 0) # Act
    Parameter('Bad_0'   , 0) # Ena
    Parameter('Bcl2_0'  , 30) # Bcl2
    Parameter('Bax_0'   , 60) # InBax

    alias_model_components()

    Initial(Bid(bf=None, state='T'), Bid_0)
    Initial(Bad(bf=None, state='M'), Bad_0)
    Initial(Bax(bf=None, s1=None, s2=None, state='C'), Bax_0)
    Initial(Bcl2(bf=None), Bcl2_0)

    # Aliases for Bax state
    inactive_monomer = {'state':'C', 's1': None, 's2': None}
    active_monomer = {'state':'A', 's1': None, 's2': None}

    # One-step "kiss-and-run" activation of Bax by tBid
    catalyze_one_step_reversible(
        Bid(state='T', bf=None), Bax(bf=None, **inactive_monomer),
        Bax(bf=None, **active_monomer), [1e-3, 1e-3])

    # Bcl2 binds tBid and Bad (a sensitizer) but not Bax
    bind_table([[                         Bcl2],
                [Bid(state='T'),  (1e-4, 1e-3)],
                [Bad(state='M'),  (1e-4, 1e-3)]])

    if pore_assembly:
        # Four Bax monomers cooperatively bind to form a tetramer
        assemble_pore_spontaneous(Bax(state='A', bf=None), [4*1e-3, 1e-3])

def cui2008_direct():
    # Build on the direct model from Chen et al. (2007) FEBS Lett. by:
    chen2007FEBS_direct(pore_assembly=False)
    alias_model_components()

    # 1. Overriding some parameter values,
    one_step_BidT_BaxC_to_BidT_BaxA_kf.value = 0.0005 # originally 0.001
    bind_BidT_Bcl2_kf.value = 0.001                   # originally 0.0001

    # 2. Adding reversible displacement reactions,
    active_monomer = {'state':'A', 's1': None, 's2': None}
    inactive_monomer = {'state':'C', 's1': None, 's2': None}
    displace_reversibly(Bax(**active_monomer), Bid(state='T'), Bcl2,
                        [0.005, 0.001])
    displace_reversibly(Bad(state='M'), Bid(state='T'), Bcl2,
                        [0.0001, 0.001])
    displace_reversibly(Bad(state='M'), Bax(active_monomer), Bcl2,
                        [0.0001, 0.005])

    # 3. Adding Bax dimerization
    active_unbound = {'state': 'A', 'bf': None}
    Rule('dimerize_Bax',
         Bax(s1=None, s2=None, **active_unbound) +
         Bax(s1=None, s2=None, **active_unbound) <>
         Bax(s1=1, s2=2, **active_unbound) % Bax(s1=2, s2=1, **active_unbound),
         Parameter('dimerize_Bax_kf', 2*0.0002),
         Parameter('dimerize_Bax_kr', 0.02))

    # 4. Adding synthesis and degradation reactions
    Bax2 = Bax(s1=1, s2=2) % Bax(s1=2, s2=1)
    synthesize_degrade_table(
       [[Bax(bf=None, **inactive_monomer),   0.06,  0.001],
        [Bax(bf=None, **active_monomer),     None,  0.001],
        [Bid(state='T', bf=None),           0.001,  0.001],
        [Bcl2(bf=None),                      0.03,  0.001],
        [Bid(state='T', bf=1) % Bcl2(bf=1),  None,  0.005],
        [Bax(bf=1) % Bcl2(bf=1),             None,  0.005],
        [Bad(state='M', bf=None),           0.001,  0.001],
        [Bad(bf=1) % Bcl2(bf=1),             None,  0.005],
        [Bax2,                               None, 0.0005]])


def cui2008_direct1(**kwargs):
    alias_model_components()

    # Build on the base "direct" model...
    cui2008_direct(**kwargs)

    # ...by adding inhibition of Bax by Bcl2
    bind(Bax(state='A', s1=None, s2=None), Bcl2, [0.005, 0.001])

def cui2008_direct2(**kwargs):
    alias_model_components()

    # Build on the "direct 1" model...
    cui2008_direct1(**kwargs)

    # By adding auto-activation of Bax
    catalyze_one_step(Bax(state='A', bf=None, s1=None, s2=None),
                      Bax(state='C', bf=None, s1=None, s2=None),
                      Bax(state='A', bf=None, s1=None, s2=None),
                      0.0002)

def howells2011():
    chen2007BiophysJ()
    # Reactions regarding Bad
    #bind(Bcl2, Bad, [1, 1])
    #displace(tBid(s=1) % Bcl2(s=1), Bad, Bad(s=1) % Bcl2(s = 1), k)
    #two_state_equilibrium(Bad(loc='M'), Bad(loc='C'), [kf, kr])
    #transition_table([[Bad(s=1) % Bcl2(s=1), Bad()????
    #                  [pBad, pBad:14-3-3, k]
    #                  [pBad:14-3-3, Bad, kBadRel],
    #                  [Bad, pBad, kBADphos1]])


