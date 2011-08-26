from pysb import *
from pysbhelperfuncs import *


Model()
# Monomers for all modules (including imported modules)
# =====================================================
Monomer('L', ['bf']) # Ligand
Monomer('R', ['bf']) # Receptor
Monomer('DISC', ['bf']) # DISC
Monomer('flip', ['bf']) # flip
Monomer('C8', ['bf', 'state'], {'state':['pro', 'A']}) # Csp 8, states: pro, active
Monomer('BAR', ['bf']) # BAR
Monomer('Bid', ['bf', 'state'], {'state':['U', 'T', 'M']}) # Bid, states: Untruncated, Truncated, truncated+Membrane
Monomer('Bax', ['bf', 'bh3', 'd2', 'state'], {'state':['C', 'M', 'A']}) # Bax, states: Cytoplasm, Mitochondria, Active
Monomer('Bak', ['bf', 'bh3', 'd2', 'state'], {'state':['M', 'A']}) # Bax, states: inactive+Membrane, Active
Monomer('Bcl2', ['bf']) # Bcl2, states: Cytoplasm, Mitochondria
Monomer('BclxL', ['bf', 'state'], {'state':['C', 'M']}) # BclxL states: cytoplasm, mitochondris
Monomer('Mcl1', ['bf']) 
Monomer('Bad', ['bf']) 
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

# EARM 1.0 Parameters and Modules 
# ===============================
from earm_1_5annbidsmacparpprms import parameter_dict as kd 
import earm_1_0modules # Must be called after the Monomers and Parameters are defined

# tBID to MOMP 
# ======================
# Bcl2, Bid, Bax migration to mitochondria
# ----------------------------------------
Rule('Bid_to_mem', Bid(bf = None, state = 'T') <> Bid(bf=None, state = 'M'), kd['BID_trans'][0],kd['BID_trans'][1])
Rule('Bax_to_mem', Bax(bf = None, state = 'C') <> Bax(bf=None, state = 'M'), kd['BAX_trans'][0], kd['BAX_trans'][1])
Rule('BclxL_to_mem', BclxL(bf = None, state = 'C') <> BclxL(bf=None, state = 'M'), kd['BCLXL_trans'][0], kd['BCLXL_trans'][1])

# Mitochondrial tBid activates Bax/Bak
# Bax/Bak form pores
# ------------------------------------
#        Bax + tBid <--> Bax:tBid --> Bax* + tBid 
#        Bak + tBid <--> Bak:tBid --> Bak* + tBid
#        Bax + Bax <--> Bax:Bax + Bax <--> Bax:Bax:Bax + Bax <--> Bax:Bax:Bax:Bax
#        Bak + Bak <--> Bak:Bak + Bak <--> Bak:Bak:Bak + Bak <--> Bak:Bak:Bak:Bak
#        Bax:Bax:Bax:Bax --> BaxPore
#        Bak:Bak:Bak:Bak --> BakPore
two_step_mod(Bid(state = 'M'), Bax(state='M'), Bax(bf = None, state = 'A'), kd['BID_BAX'])
two_step_mod(Bid(state = 'M'), Bak(state='M'), Bak(bf = None, state = 'A'), kd['BID_BAK'])
# pore_assembly(Subunit, size, rates):
pore_assembly(Bax(bf=None, state='A'), 4, kd['BAX_PORE'])
pore_assembly(Bak(bf=None, state='A'), 4, kd['BAK_PORE'])

# ------------------------------------
# MOMP Inhibition
# ------------------------------------
# Bcl2 inhibitors of Bax, Bak, and Bid
# a set of simple bind reactions:
#        Inh + Act <--> Inh:Act
# ------------------------------------
simple_bind_table([[                                            Bcl2, BclxL,  Mcl1],
                   [                                              {},    {},    {}],
                   [Bid, {'state':'M'},                         True,  False, False],
                   [Bax, {'bh3':None, 'd2':None, 'state':'A'},  True,   True, False],
                   [Bak, {'bh3':None, 'd2':None, 'state':'A'}, False,   True,  True]],
                  kd['BID_BAX_BAK_inh'], model)

# Sensitizers
# Bcl2 sensitizers bind through a simple bind resction: 
#        Inh + Act <--> Inh:Act
# This goes through the list in row-major order (as it should be)
# ---------------------------------------------------------------
simple_bind_table([[           Bcl2, BclxL,  Mcl1],
                   [             {},    {},    {}],
                   [Bad,  {},  True,  True, False],
                   [NOXA, {},  False, True,  True]],
                  kd['BCLs_sens'], model)

# Import necessary modules
# ========================
# Generate the Receptor to Bid section from the EARM 1.0 module
earm_1_0modules.rec_to_bid(model, kd)
# Generate the Pore to MOMP section from the EARM 1.0 module
earm_1_0modules.pore_to_parp(model, kd)

# Initial non-zero species
# ========================
Initial(L(bf=None), L_0)
Initial(R(bf=None), R_0)
Initial(flip(bf=None), flip_0)
Initial(C8(bf=None, state='pro'), C8_0)
Initial(BAR(bf=None), BAR_0)
Initial(Bid(bf=None, state='U'), Bid_0)
Initial(Bax(bf=None, bh3=None, d2=None, state='C'), Bax_0)
Initial(Bak(bf=None, bh3=None, d2=None, state='M'), Bax_0)
Initial(Bcl2(bf=None), Bcl2_0)
Initial(BclxL (bf=None, state='C'), BclxL_0)
Initial(Mcl1(bf=None), Mcl1_0)
Initial(Bad(bf=None), Bad_0)
Initial(NOXA(bf=None), NOXA_0)
Initial(CytoC(bf=None, state='M'), CytoC_0)
Initial(Smac(bf=None, state='M'), Smac_0)
Initial(Apaf(bf=None, state='I'), Apaf_0)
Initial(C3(bf=None, state='pro'), C3_0)
Initial(C6(bf=None, state='pro'), C6_0)
Initial(C9(bf=None), C9_0)
Initial(PARP(bf=None, state='U'), PARP_0)
Initial(XIAP(bf=None), XIAP_0)

# Observables
# ===========
# Fig 4B from Albeck observes these, normalizes and inverts them
# Observe('Bid',   Bid(bf=None, state='U'))
# Observe('PARP',  PARP(bf=None, state='U'))
# Observe('Smac',  Smac(bf=None, state='mito'))
# # This is what *should* be observed???
Observe('tBid',  Bid(state='M'))
Observe('cPARP', PARP(state='C'))
Observe('cSmac', Smac(state='A'))
Observe('cSmac_n', Smac(bf=None, state='A'))
Observe('cSmac_X', Smac(bf=1, state='A') % XIAP(bf=1))

