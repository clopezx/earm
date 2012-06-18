import pysb.macros as macros
from pysb import *

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
    return macros.catalyze(enz, 'bf', sub, 'bf', product, klist)

def bind(a, b, klist):
    return macros.bind(a, site_name, b, site_name, klist)

def bind_table(table):
    return macros.bind_table(table, site_name, site_name)

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

