from pysb import *

def declare_monomers():
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

