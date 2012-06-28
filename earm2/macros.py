import pysb.macros as macros
from pysb import *
from pysb import MonomerPattern, ComplexPattern
from pysb.util import alias_model_components

site_name = 'bf'

def declare_all_monomers():
    # Monomers for all modules (including imported modules)
    # =====================================================
    Monomer('L', ['bf']) # Ligand
    Monomer('R', ['bf']) # Receptor
    Monomer('DISC', ['bf']) # DISC
    Monomer('flip', ['bf']) # flip
    Monomer('C8', ['bf', 'state'], {'state':['pro', 'A']}) # Csp 8, states: pro, active
    Monomer('BAR', ['bf']) # BAR
    Monomer('Bid', ['bf', 'state'], {'state':['U', 'T', 'M']}) # Bid, states: Untruncated, Truncated, truncated+Membrane
    Monomer('Bax', ['bf', 's1', 's2', 'state'], {'state':['C', 'M', 'A']}) # Bax, states: Cytoplasm, Mitochondria, Active
    Monomer('Bak', ['bf', 's1', 's2', 'state'], {'state':['M', 'A']}) # Bax, states: inactive+Membrane, Active
    Monomer('Bcl2', ['bf']) # Bcl2, states: Cytoplasm, Mitochondria
    Monomer('BclxL', ['bf', 'state'], {'state':['C', 'M']}) # BclxL states: cytoplasm, mitochondris
    Monomer('Mcl1', ['bf'])
    Monomer('Bad', ['bf', 'state'], {'state':['C', 'M']})
    Monomer('NOXA', ['bf'])
    Monomer('CytoC', ['bf', 'state'], {'state':['M', 'C', 'A']})
    Monomer('Smac', ['bf', 'state'], {'state':['M', 'C', 'A']})
    Monomer('Apaf', ['bf', 'state'], {'state':['I', 'A']})
    Monomer('Apop', ['bf'])
    Monomer('C3', ['bf', 'state'], {'state':['pro', 'A', 'ub']}) # Csp 3, states: pro, active, ubiquitinated
    Monomer('C6', ['bf', 'state'], {'state':['pro', 'A']}) # Csp 6, states: pro, active
    Monomer('C9', ['bf'])
    Monomer('PARP', ['bf', 'state'], {'state':['U', 'C']}) # PARP, states: uncleaved, cleaved
    Monomer('XIAP', ['bf'])

def declare_MOMP_monomers():
    Monomer('Bid', ['bf', 'state'], {'state':['U', 'T', 'M']}) # Bid, states: Untruncated, Truncated, truncated+Membrane
    Monomer('Bax', ['bf', 's1', 's2', 'state'], {'state':['C', 'M', 'A']}) # Bax, states: Cytoplasm, Mitochondria, Active
    Monomer('Bak', ['bf', 's1', 's2', 'state'], {'state':['M', 'A']}) # Bax, states: inactive+Membrane, Active
    Monomer('Bcl2', ['bf']) # Bcl2, states: Cytoplasm, Mitochondria
    Monomer('BclxL', ['bf', 'state'], {'state':['C', 'M']}) # BclxL states: cytoplasm, mitochondris
    Monomer('Mcl1', ['bf'])
    Monomer('Bad', ['bf', 'state'], {'state':['C', 'M']})
    Monomer('NOXA', ['bf'])
    Monomer('Smac', ['bf', 'state'], {'state':['M', 'C', 'A']})

def declare_all_initial_conditions(model_type):
    if model_type not in ['indirect', 'direct', 'embedded']:
        raise ValueError("model_type must be one of 'direct', 'indirect', or 'embedded'.")

    Parameter('L_0',       3000) # 3000 Ligand correspond to 50 ng/ml SuperKiller TRAIL
    Parameter('R_0'     ,   200) # 200 TRAIL receptor 
    Parameter('flip_0'  , 1.0e2) # Flip
    Parameter('C8_0'    , 2.0e4) # procaspase-8 
    Parameter('BAR_0'   , 1.0e3) # Bifunctional apoptosis regulator
    Parameter('Bid_0'   , 4.0e4) # Bid
    Parameter('BclxL_0' , 2.0e4) # cytosolic BclxL
    Parameter('Mcl1_0'  , 2.0e4) # mitochondrial Mcl1  
    Parameter('Bad_0'   , 1.0e3) # Bad
    Parameter('NOXA_0'  , 1.0e3) # NOXA
    Parameter('CytoC_0' , 5.0e5) # cytochrome c
    Parameter('Smac_0'  , 1.0e5) # Smac    
    Parameter('Apaf_0'  , 1.0e5) # Apaf-1
    Parameter('C3_0'    , 1.0e4) # procaspase-3 (pro-C3)
    Parameter('C6_0'    , 1.0e4) # procaspase-6 (pro-C6)  
    Parameter('C9_0'    , 1.0e5) # procaspase-9 (pro-C9)
    Parameter('XIAP_0'  , 1.0e5) # X-linked inhibitor of apoptosis protein  
    Parameter('PARP_0'  , 1.0e6) # C3* substrate

    if model_type == 'indirect':
        Parameter('Bax_BclxL_0', 0.8e5), # bax + bclxl
        Parameter('Bak_Mcl1_0', 0.2e5), # bak + mcl1
        Parameter('Bax_0'   , 0) # Bax
        Parameter('Bak_0'   , 0) # Bak
    else:
        Parameter('Bcl2_0'  , 2.0e4) # cytosolic Bcl2
        Parameter('Bax_0'   , 0.8e5) # Bax
        Parameter('Bak_0'   , 0.2e5) # Bak

    alias_model_components()    

    # Initial non-zero species
    # ========================
    Initial(L(bf=None), L_0)
    Initial(R(bf=None), R_0)
    Initial(flip(bf=None), flip_0)
    Initial(C8(bf=None, state='pro'), C8_0)
    Initial(BAR(bf=None), BAR_0)
    Initial(Bid(bf=None, state='U'), Bid_0)
    Initial(Bax(bf=None, s1=None, s2=None, state='C'), Bax_0)
    Initial(BclxL (bf=None, state='C'), BclxL_0)
    Initial(Mcl1(bf=None), Mcl1_0)
    Initial(NOXA(bf=None), NOXA_0)
    Initial(CytoC(bf=None, state='M'), CytoC_0)
    Initial(Smac(bf=None, state='M'), Smac_0)
    Initial(Apaf(bf=None, state='I'), Apaf_0)
    Initial(C3(bf=None, state='pro'), C3_0)
    Initial(C6(bf=None, state='pro'), C6_0)
    Initial(C9(bf=None), C9_0)
    Initial(PARP(bf=None, state='U'), PARP_0)
    Initial(XIAP(bf=None), XIAP_0)

    if model_type == 'indirect':
        # Bak is constitutively active
        Initial(Bak(bf=None, s1=None, s2=None, state='A'), Bak_0)
        # Bad starts out at the membrane
        Initial(Bad(bf=None, state='M'), Bad_0)
        # Subpopulations of Bax and Bak are in complex with inhibitors
        Initial(Bax(bf=1, s1=None, s2=None, state='A') % BclxL(bf=1, state='M'), Bax_BclxL_0)
        Initial(Bak(bf=1, s1=None, s2=None, state='A') % Mcl1(bf=1), Bak_Mcl1_0)
    else:
        Initial(Bak(bf=None, s1=None, s2=None, state='M'), Bak_0)
        Initial(Bad(bf=None, state='C'), Bad_0)
        Initial(Bcl2(bf=None), Bcl2_0) # not used in indirect

def declare_all_observables():
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

## Aliases to pysb.macros
def catalyze(enz, sub, product, klist):
    return macros.catalyze(enz, 'bf', sub, 'bf', product, klist)

def bind(a, b, klist):
    return macros.bind(a, site_name, b, site_name, klist)

def bind_table(table):
    return macros.bind_table(table, site_name, site_name)

def assemble_pore_sequential(subunit, size, klist):
    return macros.assemble_pore_sequential(subunit, 's1', 's2', size, klist)

def pore_transport(subunit, min_size, max_size, csource, cdest, klist):
    return macros.pore_transport(subunit, 's1', 's2', 'bf', min_size, max_size,
                                csource, cdest, 'bf', klist)

## Macros for the Shen models
def assemble_pore(monomer, size, pore, klist):
    monomerp = monomer()
    kf, kr = klist
    Rule('%s_simple_pore' % monomerp.monomer.name,
         monomerp + monomerp + monomerp + monomerp <> pore,
         kf , kr)

def displace_reversibly(target, lig1, lig2, klist):
    kf, kr = klist

    r = Rule('displacement',
         target({site_name:1}) % lig1({site_name:1}) + lig2({site_name:None}) <>
         target({site_name:1}) % lig2({site_name:1}) + lig1({site_name:None}),
         Parameter('displace_kf'), Parameter('displace_kr'))

def synthesize_degrade(species, ksynth, kdeg):
    ksynth = Parameter('%s_ksynth' % species().monomer.name, ksynth)
    kdeg = Parameter('%s_kdeg' % species().monomer.name, kdeg)

    Rule('synthesize_%s' % species().monomer.name,
         None >> species(), ksynth) 
    Rule('degrade_%s' % species().monomer.name,
         species() >> None, kdeg) 

## Macros for the Albeck models
def two_step_conv(sub1, sub2, product, klist, site='bf'):
    """Automation of the Sub1 + Sub2 <> Sub1:Sub2 >> Prod two-step reaction (i.e. dimerization).
    This function assumes that there is a site named 'bf' (bind site for fxn)
    which it uses by default. Site 'bf' need not be passed when calling the function.
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

def one_step_conv(sub1, sub2, product, klist, site='bf'):
    """ Bind sub1 and sub2 to form one product: sub1 + sub2 <> product.
    """
    kf, kr = klist

    macros._verify_sites(sub1, site)
    macros._verify_sites(sub2, site)

    return macros._macro_rule('convert',
                       sub1({site: None}) + sub2({site: None}) <> product,
                       klist, ['kf', 'kr']) 

def dimerize(*args):
    pass

def bind_and_convert(*args):
    pass
