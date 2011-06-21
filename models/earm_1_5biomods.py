from pysb import *
from pysbhelperfuncs import *


Model()
# Monomers
# ======================
#('bf' site can be used with automated rule gen fxns)
Monomer('L', ['bf']) # Ligand
Monomer('R', ['bf']) # Receptor
Monomer('DISC', ['bf']) # DISC
Monomer('flip', ['bf']) # flip
Monomer('C8', ['bf', 'state'], {'state':['pro', 'A']}) # Csp 8, states: pro, active
Monomer('BAR', ['bf']) # BAR
Monomer('Bid', ['bf', 'state'], {'state':['U', 'T']}) # Bid, states: untruncated, truncated
Monomer('Bax', ['bf', 'state'], {'state':['I', 'A', 'M']}) # Bax, states: inactive, active, mitochondrial
Monomer('Bax2', ['bf'])
Monomer('Bax4', ['bf'])
Monomer('Bcl2', ['bf', 'state'], {'state':['cyto', 'mito']}) # Bcl2, states: cytoplasm, mitochondria
Monomer('MitoP', ['bf', 'state'],{'state':['U', 'A']})
Monomer('CytoC', ['bf', 'state'], {'state':['mito', 'A', 'cyto']})
Monomer('Smac', ['bf', 'state'], {'state':['mito', 'A', 'cyto']})
Monomer('Apaf', ['bf', 'state'], {'state':['I', 'A']})
Monomer('Apop', ['bf'])
Monomer('C3', ['bf', 'state'], {'state':['pro', 'A', 'ub']}) # Csp 3, states: pro, active, ubiquitinated
Monomer('C6', ['bf', 'state'], {'state':['pro', 'A']}) # Csp 6, states: pro, active
Monomer('C9', ['bf'])
Monomer('PARP', ['bf', 'state'], {'state':['U', 'C']}) # PARP, states: uncleaved, cleaved
Monomer('XIAP', ['bf'])

# EARM 1.0 Parameters and Modules 
# ===============================
import earm_1_0parms
import earm_1_0modules # Must be called after the Monomers and Parameters are defined

# tBID to MOMP 
# ======================
# Bcl2 effectors to Pore
# ----------------------
#        Bax + tBid <--> Bax:tBid --> aBax + tBid 
#        aBax <-->  MBax 
#        MBax + MBax <-->  Bax2
#        Bax2 + Bax2 <-->  Bax4
#        Bax4 + MitoP <-->  Bax4:MitoP -->  AMitoP  
# ---------------------
twostepmod(Bid(state = 'T'), Bax(state='I'), Bax(bf = None, state = 'A'), kbidbaxf, kbidbaxr, kbidbaxc)
Rule('baxCtoM', Bax(bf = None, state = 'A') <> Bax(bf=None, state = 'M'), kbaxCbaxMf, kbaxCbaxMr)
simpledim(Bax(state='M'), Bax2(bf=None), kbaxdimf, kbaxdimr)
simpledim(Bax2(), Bax4(bf = None), kbaxtetf, kbaxtetr)
twostepconv(Bax4(), MitoP(state='U'), MitoP(bf = None, state='A'), kbax4poref,kbax4porer,kbax4porec) 
# -----------------------
# Bcl2 inhibitors to Pore
# -----------------------
#        tBid + Bcl2c <-->  tBid:Bcl2c INHIBITION IN CYTO ONLY
#        MBax + Bcl2 <-->  MBax:Bcl2  
#        Bax2 + Bcl2 <-->  Bax2:Bcl2  
#        Bax4 + Bcl2 <-->  Bax4:Bcl2  
# -----------------------
simplebind(Bid(state='T'), Bcl2(state='cyto'), kbidbcl2f, kbidbcl2r)
simplebind(Bax(state='M'), Bcl2(state='mito'), kbaxMbcl2Mf, kbaxMbcl2Mr)
simplebind(Bax2(), Bcl2(state='mito'), kbax2Mbcl2Mf, kbax2Mbcl2Mr)
simplebind(Bax4(), Bcl2(state='mito'), kbax4Mbcl2Mf, kbax4Mbcl2Mr)

# Import necessary modules
# ========================
# Generate the Receptor to Bid section from the EARM 1.0 module
earm_1_0modules.rec_to_bid(model)
# Generate the Pore to MOMP section from the EARM 1.0 module
earm_1_0modules.pore_to_parp(model)

# Initial amounts
# ===============
# including species from modules!
# -------------------------------
Initial(L(bf=None), L_0)
Initial(R(bf=None), R_0)
Initial(flip(bf=None), flip_0)
Initial(C8(bf=None, state='pro'), C8_0)
Initial(BAR(bf=None), BAR_0)
Initial(C3(bf=None, state='pro'), C3_0)
Initial(C6(bf=None, state='pro'), C6_0)
Initial(XIAP(bf=None), XIAP_0)
Initial(PARP(bf=None, state='U'), PARP_0)
Initial(Bid(bf=None, state='U'), Bid_0)
Initial(Bcl2(bf=None, state='cyto'), Bcl2_cyto_0)
Initial(Bax(bf=None, state='I'), Bax_0)
Initial(Bcl2(bf=None, state='mito'), Bcl2_mito_0)
Initial(MitoP(bf=None, state='U'), MitoP_0)
Initial(CytoC(bf=None, state='mito'), CytoC_0)
Initial(Smac(bf=None, state='mito'), Smac_0)
Initial(C9(bf=None), C9_0)
Initial(Apaf(bf=None, state='I'), Apaf_0)

# Observables
# ===========
# Fig 4B from Albeck observes these, normalizes and inverts them
Observe('Bid',   Bid(bf=None, state='U'))
Observe('PARP',  PARP(bf=None, state='U'))
Observe('Smac',  Smac(bf=None, state='mito'))
# This is what *should* be observed???
Observe('tBid',  Bid(state='T'))
Observe('cPARP', PARP(state='C'))
Observe('cSmac', Smac(state='cyto'))

