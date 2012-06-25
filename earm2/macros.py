from pysb import *
import pysb.macros
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
    Monomer('tBid', ['bf'])
    Monomer('Bax', ['bf', 'state'], {'state':['C', 'A']})
    Monomer('Bcl2', ['bf'])
    Monomer('Bad', ['bf'])
    Monomer('Pore')

def declare_all_initial_conditions():
    alias_model_components()
    # TODO: make local macro "initial" which takes a float and automatically create a Parameter
    Initial(L(bf=None), L_0)
    Initial(R(bf=None), R_0)
    Initial(flip(bf=None), flip_0)
    Initial(C8(bf=None, state='pro'), C8_0)
    Initial(BAR(bf=None), BAR_0)
    Initial(Bid(bf=None, state='U'), Bid_0)
    Initial(Bax(bf=None, s1=None, s2=None, state='C'), Bax_0)
    Initial(Bak(bf=None, s1=None, s2=None, state='M'), Bak_0)
    Initial(Bcl2(bf=None), Bcl2_0)
    Initial(BclxL (bf=None, state='C'), BclxL_0)
    Initial(Mcl1(bf=None), Mcl1_0)
    Initial(Bad(bf=None, state='C'), Bad_0)
    Initial(NOXA(bf=None), NOXA_0)
    Initial(CytoC(bf=None, state='M'), CytoC_0)
    Initial(Smac(bf=None, state='M'), Smac_0)
    Initial(Apaf(bf=None, state='I'), Apaf_0)
    Initial(C3(bf=None, state='pro'), C3_0)
    Initial(C6(bf=None, state='pro'), C6_0)
    Initial(C9(bf=None), C9_0)
    Initial(PARP(bf=None, state='U'), PARP_0)
    Initial(XIAP(bf=None), XIAP_0)    

def catalyze_one_step_reversible(enz, sub, prod, klist):
    kcat, krev = klist

    Rule('catalyze1',
         enz({site_name:None}) + sub({site_name:None}) >>
         enz({site_name:None}) + prod({site_name:None}),
         Parameter('catalyze1_kcat', kcat))

    Rule('catalyze1_rev',
         prod({site_name:None}) >> sub({site_name:None}),
         Parameter('cat1_rev_krev', krev))

def catalyze(enz, sub, product, klist):
    return pysb.macros.catalyze(enz, 'bf', sub, 'bf', product, klist)

def bind(a, b, klist):
    return pysb.macros.bind(a, site_name, b, site_name, klist)

def bind_table(table):
    return pysb.macros.bind_table(table, site_name, site_name)

def assemble_pore(monomer, size, pore, klist):
    monomerp = monomer()
    kf, kr = klist
    Rule('%s_simple_pore' % monomerp.monomer.name,
         monomerp + monomerp + monomerp + monomerp <> pore,
         kf , kr)

def assemble_pore_sequential(monomer, size, klist):
    return pysb.macros.assemble_pore_sequential(monomer, 's1', 's2', size,
                                                klist)

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

def translocate(sub, state1, state2, kf=1e-2, kr=1e-2):
    return pysb.macros.two_state_equilibrium(sub(state=state1),
                                             sub(state=state2), [kf, kr])
