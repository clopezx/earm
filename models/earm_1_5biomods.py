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
Monomer('Bid', ['bf', 'state'], {'state':['U', 'T', 'M']}) # Bid, states: Untruncated, Truncated, truncated+Membrane
Monomer('Bax', ['bf', 'state'], {'state':['C', 'M', 'A']}) # Bax, states: Cytoplasm, Mitochondria, Active
Monomer('Bak', ['bf', 'state'], {'state':['M', 'A']}) # Bax, states: inactive+Membrane, Active
Monomer('Bcl2', ['bf', 'state'], {'state':['C', 'M']}) # Bcl2, states: Cytoplasm, Mitochondria
Monomer('Mcl1', ['bf']) 
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
import earm_1_5parms
import earm_1_0modules # Must be called after the Monomers and Parameters are defined

# tBID to MOMP 
# ======================
# Bcl2, Bid, Bax migration to mitochondria
# ----------------------------------------
Rule('Bax_to_mem', Bax(bf = None, state = 'C') <> Bax(bf=None, state = 'M'), kbaxCbaxMf, kbaxCbaxMr)
Rule('Bcl2_to_mem', Bcl2(bf = None, state = 'C') <> Bcl2(bf=None, state = 'M'), kbaxCbaxMf, kbaxCbaxMr)
Rule('Bid_to_mem', Bid(bf = None, state = 'T') <> Bid(bf=None, state = 'M'), kbaxCbaxMf, kbaxCbaxMr)
# Mitochondrial tBid activates Bax/Bak
# Bax/Bak form pores
# ------------------------------------
#        Bax + tBid <--> Bax:tBid --> Bax* + tBid 
#        Bak + tBid <--> Bak:tBid --> Bak* + tBid
#        Bax + Bax <--> Bax:Bax + Bax <--> Bax:Bax:Bax + Bax <--> Bax:Bax:Bax:Bax
#        Bak + Bak <--> Bak:Bak + Bak <--> Bak:Bak:Bak + Bak <--> Bak:Bak:Bak:Bak
#        Bax:Bax:Bax:Bax --> BaxPore
#        Bak:Bak:Bak:Bak --> BakPore
twostepmod(Bid(state = 'M'), Bax(state='M'), Bax(bf = None, state = 'A'), kbidbaxf, kbidbaxr, kbidbaxc)
twostepmod(Bid(state = 'M'), Bak(state='M'), Bax(bf = None, state = 'A'), kbidbaxf, kbidbaxr, kbidbaxc)
oligomerize(Bax(state='A'), 4) # oligomerize to tetramer
oligomerize(Bak(state='A'), 4) # oligomerize to tetramer
Rule('tetraBax_to_pore', 
     Bax(b1=None, b2=1)%Bax(b1=1, b2=2)%Bax(b1=2,b2=3)%Bax(b1=3,b2=None) <> BaxPore(), kf, kr)
Rule('tetraBak_to_pore', 
     Bak(b1=None, b2=1)%Bak(b1=1, b2=2)%Bak(b1=2,b2=3)%Bak(b1=3,b2=None) <> BakPore(), kf, kr)
# ------------------------------------
# MOMP Inhibition
# ------------------------------------
#        aBax <-->  MBax 
#        MBax + MBax <-->  Bax2
#        Bax2 + Bax2 <-->  Bax4
#        Bax4 + MitoP <-->  Bax4:MitoP -->  AMitoP  
# ---------------------
simplebind(Bid(state='M'), Bcl2(state='M'), kbidbcl2f, kbidbcl2r)
simplebind(Bax(), Bcl2(state='M'), kbaxMbcl2Mf, kbaxMbcl2Mr) # CHECK: only membrane-bound heterodimers are formed in the rulz
simplebind(Bak(), Mcl1(state='M'), kbaxMbcl2Mf, kbaxMbcl2Mr) # CHECK: only membrane-bound heterodimers are formed in the rulz

XXsimplebind(Bax2(), Bcl2(state='mito'), kbax2Mbcl2Mf, kbax2Mbcl2Mr)
XXsimplebind(Bax4(), Bcl2(state='mito'), kbax4Mbcl2Mf, kbax4Mbcl2Mr)

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

