from pysb import *
from macros import *
from pysb.macros import catalyze_one_step_reversible, catalyze_one_step, \
                        synthesize_degrade_table, degrade, equilibrate
from pysb.util import alias_model_components

transloc_rates = [1e-2, 1e-2]

# Useful aliases for Bax state
inactive_monomer = {'state':'C', 's1': None, 's2': None}
active_monomer = {'state':'A', 's1': None, 's2': None}

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

def chen2007BiophysJ(pore_assembly=True):
    # TODO: change all initial conditions and param values to Molar
    Parameter('Bid_0'   , 0) # Bid
    Parameter('Bcl2_0'  , 1e-1) # Mitochondrial Bcl2
    Parameter('Bax_0'   , 2e-1) # Bax

    alias_model_components()

    Initial(Bid(bf=None, state='T'), Bid_0)
    Initial(Bax(bf=None, s1=None, s2=None, state='C'), Bax_0)
    Initial(Bcl2(bf=None), Bcl2_0)

    # Aliases for Bax state

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

    if pore_assembly:
        # Four Bax monomers cooperatively bind to form a tetramer
        assemble_pore_spontaneous(Bax(state='A', bf=None), [2, 0])

def chen2007FEBS_indirect(pore_assembly=True):
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
    bind_table([[                                Bcl2],
                [Bid(state='T'),         (1e-4, 1e-3)],
                [Bax(inactive_monomer),  (1e-4, 1e-3)]])

    if pore_assembly:
        # Four "inactive" Bax monomers cooperatively bind to form a tetramer
        assemble_pore_spontaneous(Bax(state='C', bf=None), [4*1e-3, 1e-3])

def chen2007FEBS_direct(pore_assembly=True):
    # TODO: change all initial conditions and param values to Molar
    # Initial conditions
    Parameter('Bid_0'   , 0) # Act
    Parameter('Bad_0'   , 0) # Ena
    Parameter('Bcl2_0'  , 30) # Bcl2
    Parameter('Bax_0'   , 60) # InBax
    alias_model_components()
    Initial(Bid(bf=None, state='T'), Bid_0)
    Initial(Bad(bf=None, state='M', serine='U'), Bad_0)
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

    # 2. Adding a Bad-for-Bid displacement reaction,
    displace_reversibly(Bad(state='M'), Bid(state='T'), Bcl2,
                        [0.0001, 0.001])

    # 3. Adding Bax dimerization,
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
       [[Bax(bf=None, **inactive_monomer),      0.06,  0.001],
        [Bax(bf=None, **active_monomer),        None,  0.001],
        [Bid(state='T', bf=None),              0.001,  0.001],
        [Bcl2(bf=None),                         0.03,  0.001],
        [Bid(state='T', bf=1) % Bcl2(bf=1),     None,  0.005],
        [Bad(state='M', bf=None, serine='U'),  0.001,  0.001],
        [Bad(bf=1) % Bcl2(bf=1),                None,  0.005],
        [Bax2,                                  None, 0.0005]])

def cui2008_direct1():
    alias_model_components()

    # Build on the base "direct" model...
    cui2008_direct()

    # ...by adding inhibition of Bax by Bcl2,
    bind(Bax(state='A', s1=None, s2=None), Bcl2, [0.005, 0.001])

    # ...associated displacement reactions
    displace_reversibly(Bax(active_monomer), Bid(state='T'), Bcl2,
                        [0.005, 0.001])
    displace_reversibly(Bad(state='M'), Bax(active_monomer), Bcl2,
                        [0.0001, 0.005])

    # ...and degradation of the active Bax:Bcl2 complex
    degrade(Bax(bf=1) % Bcl2(bf=1), 0.005)

def cui2008_direct2():
    alias_model_components()

    # Build on the "direct 1" model...
    cui2008_direct1()

    # By adding simultaneous auto-activation and dimerization of Bax
    Rule('Bax_autoactivation_dimerization',
        Bax(state='A', bf=None, s1=None, s2=None) +
        Bax(state='C', bf=None, s1=None, s2=None) >>
        Bax(state='A', bf=None, s1=1, s2=2) %
        Bax(state='A', bf=None, s1=2, s2=1),
        Parameter('Bax_autoactivation_dimerization_k', 0.0002))

# TODO: All parameter values
def howells2011():
    # Build on the model from Chen et al. (2007) Biophys J:
    chen2007BiophysJ()

    # Add initial condition for Bad
    Parameter('Bad_0', 0) # Ena
    alias_model_components()
    Initial(Bad(bf=None, state='M', serine='U'), Bad_0)

    # Translocation equilibrium between unphosphorylated cytosolic and
    # mitochondrial Bad
    equilibrate(Bad(state='C', serine='U', bf=None),
                Bad(state='M', serine='U', bf=None), [1, 1])

    # Bad binds Bcl2
    bind(Bad(state='M'), Bcl2, [1, 1])

    # Bad displaces tBid from Bcl2
    displace(Bad(state='M'), Bid(state='T'), Bcl2, 1)

    # Phosphorylation of Bad
    Rule('phosphorylate_BadCU_to_BadCP',     # Cytosolic Bad
         Bad(state='C', serine='U') >> Bad(state='C', serine='P'),
         Parameter('phosphorylate_BadC_k', 1))
    Rule('phosphorylate_BadMU_to_BadCP',     # Mitochondrial Bad
         Bad(state='M', serine='U', bf=None) >>
         Bad(state='C', serine='P', bf=None),
         Parameter('phosphorylate_BadM_k', 1))
    Rule('phosphorylate_BadMUBcl2_to_BadCP', # Mitochondrial Bad:Bcl2
         Bad(state='M', serine='U', bf=1) % Bcl2(bf=1) >>
         Bad(state='C', serine='P', bf=None) + Bcl2(bf=None),
         Parameter('phosphorylate_BadMUBcl2_to_BadCP_k', 1))

    # Sequester phospho-Bad by "binding" 14-3-3 domains
    Rule('sequester_BadCP_to_BadC1433',
         Bad(state='C', serine='P') >> Bad(state='C', serine='B'),
         Parameter('sequester_BadCP_to_BadC1433_k', 1))

    # Release of Bad from 14-3-3 domains
    Rule('release_BadC1433_to_BadCU',
         Bad(state='C', serine='B') >> Bad(state='C', serine='U'),
         Parameter('release_BadC1433_to_BadCU_k', 1))
