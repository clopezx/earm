# EARM 2.0 INDIRECT
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
Monomer('Bax', ['bf', 's1', 's2', 'state'], {'state':['C', 'M', 'A']}) # Bax, states: Cytoplasm, Mitochondria, Active
Monomer('Bak', ['bf', 's1', 's2', 'state'], {'state':['M', 'A']}) # Bax, states: inactive+Membrane, Active
Monomer('BclxL', ['bf', 'state'], {'state':['C', 'M']}) # BclxL states: cytoplasm, mitochondris
Monomer('Bcl2', ['bf']) 
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
Monomer('Bid', ['bf', 'state'], {'state':['U', 'T', 'M']}) # NOT USED. Kept here simply to use the trail_to_bid module

# EARM 1.0 Parameters and Modules 
# ===============================
from earm_indirect_mem_parms import parameter_dict as kd 
import earm_1_0modules # Must be called after the Monomers and Parameters are defined

# tBID to MOMP 
# ======================
# Bcl2, Bid, Bax migration to mitochondria
# ----------------------------------------
free_Bax = Bax(bf=None, s1=None, s2=None)
Rule('Bax_to_mem', free_Bax(state = 'C') <> free_Bax(state = 'A'), *kd['BAX_trans'])
Rule('BclxL_to_mem', BclxL(bf = None, state = 'C') <> BclxL(bf = None, state = 'M'), *kd['BCLXL_trans'])

# Bax/Bak form pores
# ------------------------------------
# ringp_assembly(Subunit, size, rates):
ringp_assembly(Bax(bf=None, state='A'), 4, kd['BAX_PORE'])
ringp_assembly(Bak(bf=None, state='A'), 4, kd['BAK_PORE'])

# ------------------------------------
# MOMP Inhibition
# ------------------------------------
# Bcl2 inhibitors of Bax, Bak, and Bid
# ------------------------------------
simple_bind_table([[                                                   BclxL,  Mcl1,  Bcl2],
                   [                                           {'state':'M'},    {},    {}], 
                   [Bid, {'state':'T'},                                 True,  True,  True], #NOTE: MOVE TO SENSITIZER 
                   [Bax, {'s1':None, 's2':None, 'state':'A'},           True, False,  True],
                   [Bak, {'s1':None, 's2':None, 'state':'A'},           True,  True, False]],
                  kd['BID_BAX_BAK_inh'], model)

# Bcl2 sensitizers 
# ---------------------------------------------------------------
simple_bind_table([[                  BclxL,  Mcl1,  Bcl2], #Note: according to Andrews, these should all be at mitochondria
                   [          {'state':'M'},    {},    {}],
                   [Bad,  {},          True, False,  True],
                   [NOXA, {},         False,  True, False]],
                  kd['BCLs_sens'], model)

# CytC, Smac release
# ----------------------
# ringp_transport(Subunit, Source, Dest, min_size, max_size, rates):
ringp_transport(Bax(bf=None), CytoC(state='M'), CytoC(state='C'), 4, 4, kd['BAX_CYTC']) 
ringp_transport(Bax(bf=None),  Smac(state='M'),  Smac(state='C'), 4, 4, kd['BAX_SMAC']) 
ringp_transport(Bak(bf=None), CytoC(state='M'), CytoC(state='C'), 4, 4, kd['BAK_CYTC'])
ringp_transport(Bak(bf=None),  Smac(state='M'),  Smac(state='C'), 4, 4, kd['BAK_SMAC'])

# --------------------------------------
# CytC and Smac activation after release
# --------------------------------------
Rule('act_cSmac',  Smac(bf=None, state='C') <> Smac(bf=None, state='A'), kd['SMAC_ACT'][0], kd['SMAC_ACT'][1])
Rule('act_cCytoC', CytoC(bf=None, state='C') <> CytoC(bf=None, state='A'), kd['CYTOC_ACT'][0], kd['CYTOC_ACT'][1])

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
Initial(Bax(bf=None, s1=None, s2=None, state='C'), Bax_0)
Initial(Bax(bf=1, s1=None, s2=None, state='A') % Bcl2(bf=1), Bax_Bcl2_0)
Initial(Bak(bf=None, s1=None, s2=None, state='A'), Bak_0)
Initial(Bak(bf=1, s1=None, s2=None, state='A') % Mcl1(bf=1), Bak_Mcl1_0)
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
Observable('tBid',  Bid(state='T'))  # 1 no state M in this case!
Observable('aSmac', Smac(state='A')) # 2
Observable('cPARP', PARP(state='C')) # 3
#Observe('BaxBclxl',  Bax(bf=1, s1=None, s2=None, state='A') % BclxL(bf=1, state='M'))
#Observe('BakMcl1', Bak(bf=1, s1=None, s2=None, state='A') % Mcl1(bf=1))
