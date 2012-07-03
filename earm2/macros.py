import pysb.macros as macros
from pysb import *
from pysb import MonomerPattern, ComplexPattern
from pysb.util import alias_model_components

# The default site name to be used for binding reactions
site_name = 'bf'

## Monomer declarations ========================
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

    ligand_to_c8_monomers()
    momp_monomers()
    apaf1_to_parp_monomers()

def ligand_to_c8_monomers():
    """ Declares ligand, receptor, DISC, Flip, Bar and Caspase 8.

    The package variable site_name specifies the name of the site to be used
    for all binding reactions.

    The 'state' site denotes various localization and/or activity states of a
    Monomer, with 'C' denoting cytoplasmic localization and 'M' mitochondrial
    localization.
    """

    Monomer('L', [site_name]) # Ligand
    Monomer('R', [site_name]) # Receptor
    Monomer('DISC', [site_name]) # Death-Inducing Signaling Complex
    Monomer('flip', [site_name])
    # Caspase 8, states: pro, Active
    Monomer('C8', [site_name, 'state'], {'state':['pro', 'A']})
    Monomer('BAR', [site_name])

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

def apaf1_to_parp_monomers():
    """ Declares CytochromeC, Smac, Apaf-1, the Apoptosome, Caspases 3, 6, 9,
    XIAP and PARP.

    The package variable site_name specifies the name of the site to be used
    for all binding reactions.

    The 'state' site denotes various localization and/or activity states of a
    Monomer, with 'C' denoting cytoplasmic localization and 'M' mitochondrial
    localization.
    """

    # Cytochrome C
    Monomer('Apaf', [site_name, 'state'], {'state':['I', 'A']}) # Apaf-1
    Monomer('Apop', [site_name]) # Apoptosome (activated Apaf-1 + caspase 9)
    # Csp 3, states: pro, active, ubiquitinated
    Monomer('C3', [site_name, 'state'], {'state':['pro', 'A', 'ub']})
    # Caspase 6, states: pro-, Active
    Monomer('C6', [site_name, 'state'], {'state':['pro', 'A']})
    Monomer('C9', [site_name]) # Caspase 9
    # PARP, states: Uncleaved, Cleaved
    Monomer('PARP', [site_name, 'state'], {'state':['U', 'C']})
    Monomer('XIAP', [site_name]) # X-linked Inhibitor of Apoptosis Protein

## Initial condition declarations =============
def all_initial_conditions():
    """Declare initial conditions for the full extrinsic apoptosis model.
    """

    ligand_to_c8_initial_conditions()
    momp_initial_conditions(bid_state='U')
    apaf1_to_parp_initial_conditions()

def ligand_to_c8_initial_conditions():
    """Declare initial conditions for ligand, receptor, Flip, C8, and Bar.
    """

    Parameter('L_0',       3000) # 3000 Ligand corresponds to 50 ng/ml SK-TRAIL
    Parameter('R_0'     ,   200) # 200 TRAIL receptor
    Parameter('flip_0'  , 1.0e2) # Flip
    Parameter('C8_0'    , 2.0e4) # procaspase-8
    Parameter('BAR_0'   , 1.0e3) # Bifunctional apoptosis regulator

    alias_model_components()

    Initial(L(bf=None), L_0)
    Initial(R(bf=None), R_0)
    Initial(flip(bf=None), flip_0)
    Initial(C8(bf=None, state='pro'), C8_0)
    Initial(BAR(bf=None), BAR_0)

def momp_initial_conditions(bid_state='U'):
    """Declare initial conditions for Bcl-2 family proteins, Cyto c, and Smac.

    Parameters
    ==========

    """
    Parameter('Bid_0'   , 4.0e4) # Bid
    Parameter('BclxL_0' , 2.0e4) # cytosolic BclxL
    Parameter('Mcl1_0'  , 2.0e4) # Mitochondrial Mcl1
    Parameter('Bcl2_0'  , 2.0e4) # Mitochondrial Bcl2
    Parameter('Bad_0'   , 1.0e3) # Bad
    Parameter('NOXA_0'  , 1.0e3) # NOXA
    Parameter('CytoC_0' , 5.0e5) # cytochrome c
    Parameter('Smac_0'  , 1.0e5) # Smac
    Parameter('Bax_0'   , 0.8e5) # Bax
    Parameter('Bak_0'   , 0.2e5) # Bak

    alias_model_components()

    Initial(Bid(bf=None, state=bid_state), Bid_0)
    Initial(Bad(bf=None, state='C'), Bad_0)
    Initial(Bax(bf=None, s1=None, s2=None, state='C'), Bax_0)
    Initial(Bak(bf=None, s1=None, s2=None, state='M'), Bak_0)
    Initial(Bcl2(bf=None), Bcl2_0)
    Initial(BclxL (bf=None, state='C'), BclxL_0)
    Initial(Mcl1(bf=None, state='M'), Mcl1_0)
    Initial(NOXA(bf=None, state='C'), NOXA_0)
    Initial(CytoC(bf=None, state='M'), CytoC_0)
    Initial(Smac(bf=None, state='M'), Smac_0)

def apaf1_to_parp_initial_conditions():
    """Declare initial conditions for CytoC, Smac, Apaf-1, Apoptosome, caspases
       3, 6, and 9, XIAP, and PARP.
    """

    Parameter('Apaf_0'  , 1.0e5) # Apaf-1
    Parameter('C3_0'    , 1.0e4) # procaspase-3 (pro-C3)
    Parameter('C6_0'    , 1.0e4) # procaspase-6 (pro-C6)
    Parameter('C9_0'    , 1.0e5) # procaspase-9 (pro-C9)
    Parameter('XIAP_0'  , 1.0e5) # X-linked inhibitor of apoptosis protein
    Parameter('PARP_0'  , 1.0e6) # C3* substrate

    alias_model_components()

    Initial(Apaf(bf=None, state='I'), Apaf_0)
    Initial(C3(bf=None, state='pro'), C3_0)
    Initial(C6(bf=None, state='pro'), C6_0)
    Initial(C9(bf=None), C9_0)
    Initial(PARP(bf=None, state='U'), PARP_0)
    Initial(XIAP(bf=None), XIAP_0)

## Observables declarations ===================
def all_observables():
    alias_model_components()
    # Observables
    # ===========
    # Fig 4B from Albeck observes these, normalizes and inverts them
    # Observe('Bid',   Bid(bf=None, state='U'))
    # Observe('PARP',  PARP(bf=None, state='U'))
    # Observe('Smac',  Smac(bf=None, state='mito'))
    # # This is what *should* be observed???
    Observable('mBid',  Bid(state='M'))
    Observable('cSmac', Smac(state='A'))
    Observable('cPARP', PARP(state='C'))

## Aliases to pysb.macros =====================
# TODO use site_name for all of these
def catalyze(enz, sub, product, klist):
    return macros.catalyze(enz, site_name, sub, site_name, product, klist)

def bind(a, b, klist):
    return macros.bind(a, site_name, b, site_name, klist)

def bind_table(table):
    return macros.bind_table(table, site_name, site_name)

def assemble_pore_sequential(subunit, size, klist):
    return macros.assemble_pore_sequential(subunit, 's1', 's2', size, klist)

def pore_transport(subunit, min_size, max_size, csource, cdest, klist):
    return macros.pore_transport(subunit, 's1', 's2', 'bf', min_size, max_size,
                                csource, 'bf', cdest, klist)

## Macros for the Shen models
def assemble_pore_spontaneous(subunit, klist):

    def pore_rule_name(rule_expression):
        react_p = rule_expression.reactant_pattern
        mp = react_p.complex_patterns[0].monomer_patterns[0]
        subunit_name = macros._monomer_pattern_label(mp)
        pore_name = mp.monomer.name
        return '%s_to_%s%d' % (subunit_name, mp.monomer.name, 4)

    free_subunit = subunit(s1=None, s2=None)
    macros._macro_rule('spontaneous_pore',
        free_subunit + free_subunit + free_subunit + free_subunit <>
        subunit(s1=1, s2=4) % subunit(s1=2, s2=1) % \
        subunit(s1=3, s2=2) % subunit(s1=4, s2=3),
        klist, ['kf', 'kr'], name_func=pore_rule_name)

def displace_reversibly(target, lig1, lig2, klist):
    """Generate displacement reaction T:L1 + L2 <> L1 + T:L2
    """
    return macros._macro_rule('displace',
         target({site_name:1}) % lig1({site_name:1}) + lig2({site_name:None}) <>
         target({site_name:1}) % lig2({site_name:1}) + lig1({site_name:None}),
         klist, ['kf', 'kr'])

def synthesize_degrade(species, ksynth, kdeg):
    ksynth = Parameter('%s_ksynth' % species().monomer.name, ksynth)
    kdeg = Parameter('%s_kdeg' % species().monomer.name, kdeg)

    Rule('synthesize_%s' % species().monomer.name,
         None >> species(), ksynth) 
    Rule('degrade_%s' % species().monomer.name,
         species() >> None, kdeg) 

## Macros for the Albeck models
def two_step_conv(sub1, sub2, product, klist, site=site_name):
    """Automation of the Sub1 + Sub2 <> Sub1:Sub2 >> Prod two-step reaction.

    Because product is created by the function, it must be fully specified.
    """

    macros._verify_sites(sub1, site)
    macros._verify_sites(sub2, site)

    components = macros._macro_rule('complex',
                             sub1({site: None}) + sub2({site: None}) <>
                             sub1({site: 1}) % sub2({site: 1}),
                             klist[0:2], ['kf', 'kr'])
    components |= macros._macro_rule('convert',
                              sub1({site: 1}) % sub2({site: 1}) >> product,
                              [klist[2]], ['kc'])
    return components

def one_step_conv(sub1, sub2, product, klist, site=site_name):
    """ Bind sub1 and sub2 to form one product: sub1 + sub2 <> product.
    """
    kf, kr = klist

    macros._verify_sites(sub1, site)
    macros._verify_sites(sub2, site)

    return macros._macro_rule('convert',
                       sub1({site: None}) + sub2({site: None}) <> product,
                       klist, ['kf', 'kr']) 
