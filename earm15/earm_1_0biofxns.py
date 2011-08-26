from pysb import *
from pysbhelperfuncs import *

# Albeck JG, Burke JM, Spencer SL, Lauffenburger DA, Sorger PK, 2008
# Modeling a Snap-Action, Variable-Delay Switch Controlling Extrinsic
# Cell Death. PLoS Biol 6(12): e299. doi:10.1371/journal.pbio.0060299
#
# http://www.plosbiology.org/article/info:doi/10.1371/journal.pbio.0060299
#
#

Model()

# Monomers
#('bf' site can be used with automated rule gen fxns)
Monomer('L', ['bf']) # Ligand
Monomer('R', ['bf']) # Receptor
Monomer('DISC', ['bf']) # DISC
Monomer('flip', ['bf']) # flip
Monomer('C8', ['bf', 'state'], {'state':['pro', 'A']}) # Csp 8, states: pro, active
Monomer('BAR', ['bf']) # BAR
Monomer('C3', ['bf', 'state'], {'state':['pro', 'A', 'ub']}) # Csp 3, states: pro, active, ubiquitinated
Monomer('C6', ['bf', 'state'], {'state':['pro', 'A']}) # Csp 6, states: pro, active
Monomer('PARP', ['bf', 'state'], {'state':['U', 'C']}) # PARP, states: uncleaved, cleaved
Monomer('Bid', ['bf', 'state'], {'state':['U', 'T']}) # Bid, states: untruncated, truncated
Monomer('Bcl2', ['bf', 'state'], {'state':['cyto', 'mito']}) # Bcl2, states: cytoplasm, mitochondria
Monomer('Bax', ['bf', 'state'], {'state':['I', 'A', 'M']}) # Bax, states: inactive, active, mitochondrial
Monomer('Bax2', ['bf'])
Monomer('Bax4', ['bf'])
Monomer('MitoP', ['bf', 'state'],{'state':['U', 'A']})
Monomer('CytoC', ['bf', 'state'], {'state':['mito', 'A', 'cyto']})
Monomer('Smac', ['bf', 'state'], {'state':['mito', 'A', 'cyto']})
Monomer('C9', ['bf'])
Monomer('Apaf', ['bf', 'state'], {'state':['I', 'A']})
Monomer('Apop', ['bf'])
Monomer('XIAP', ['bf'])

# get the parameters
import earm_1_0parms

# RECEPTOR TO tBID
# =====================
# tBID Activation Rules
# ---------------------
#        L + R <--> L:R --> DISC
#        pC8 + DISC <--> DISC:pC8 --> C8 + DISC
#        Bid + C8 <--> Bid:C8 --> tBid + C8
# ---------------------
twostepconv(L(), R(), DISC(bf = None ), klrf, klrr, klrc)
twostepmod(DISC(), C8(state='pro'), C8(bf = None, state='A'), kdiscc8f, kdiscc8r, kdiscc8c)
twostepmod(C8(state='A'), Bid(state='U'), Bid(bf = None, state='T'), kc8bidf, kc8bidr, kc8bidc)
# ---------------------
# Inhibition Rules
# ---------------------
#        flip + DISC <-->  flip:DISC  
#        C8 + BAR <--> BAR:C8 CSPS
# ---------------------
simplebind(DISC(), flip(), kflipdiscf, kflipdiscr)
simplebind(BAR(), C8(state='A'), kbarc8f, kbarc8r)
# ---------------------

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

# FXR CASPASES CLEAVE PARP 
# ========================
# CytC, Smac release
# ----------------------
#        AMito + mCytoC <-->  AMito:mCytoC --> AMito + ACytoC  
#        AMito + mSmac <-->  AMito:mSmac --> AMito + ASmac  
#        ACytoC <-->  cCytoC
#        ASmac <-->  cSmac
# ----------------------
twostepmod(MitoP(state='A'), CytoC(state='mito'), CytoC(bf = None, state='A'), kmitopcytocMf, kmitopcytocMr, kmitopcytocMc)
twostepmod(MitoP(state='A'), Smac(state='mito'), Smac(bf = None, state='A'), kmitopsmacMf, kmitopsmacMr, kmitopsmacMc)
Rule('cytocCtoM', CytoC(bf = None, state = 'A') <> CytoC(bf=None, state = 'cyto'), kcytocMcytocCf, kcytocMcytocCr)
Rule('SmacMtoSmacC', Smac(bf = None, state='A') <> Smac(bf = None, state='cyto'), ksmacMsmacCf, ksmacMsmacCr)
# ---------------------------
# Apoptosome formation
# ---------------------------
#        Apaf + cCytoC <-->  Apaf:cCytoC --> aApaf + cCytoC
#        aApaf + pC9 <-->  Apop
#        Apop + pC3 <-->  Apop:pC3 --> Apop + C3
# ---------------------------
twostepmod(CytoC(state='cyto'), Apaf(state='I'), Apaf(bf = None, state = 'A'), kcytocCapaff, kcytocCapafr, kcytocCapafc)
onestepconv(Apaf(state='A'), C9(), Apop(bf=None), kapafc9f, kapafc9r)
twostepmod(Apop(), C3(state='pro'), C3(bf = None, state='A'), kapopc3f, kapopc3r, kapopc3c)
# -----------------------------
# Apoptosome related inhibitors
# -----------------------------
#        Apop + XIAP <-->  Apop:XIAP  
#        cSmac + XIAP <-->  cSmac:XIAP  
simplebind(Apop(), XIAP(), kapopxiapf, kapopxiapr)
simplebind(Smac(state='cyto'), XIAP(), ksmacxiapf, ksmacxiapr)
# ---------------------------
# Caspase reactions (effectors, inhibitors, and loopback initiators)
# ---------------------------
#        pC3 + C8 <--> pC3:C8 --> C3 + C8 CSPS
#        pC6 + C3 <--> pC6:C3 --> C6 + C3 CSPS
#        XIAP + C3 <--> XIAP:C3 --> XIAP + C3_U CSPS
#        PARP + C3 <--> PARP:C3 --> CPARP + C3 CSPS
#        pC8 + C6 <--> pC8:C6 --> C8 + C6 CSPS
# ---------------------------
twostepmod(C8(state='A'), C3(state='pro'), C3(bf = None, state='A'), kc8c3f, kc8c3r, kc8c3c)
twostepmod(C3(state='A'), C6(state='pro'), C6(bf = None, state='A'), kc3c6f, kc3c6r, kc3c6c)
twostepmod(XIAP(), C3(state = 'A'), C3(bf = None, state = 'ub'), kxiapc3f, kxiapc3r, kxiapc3c)
twostepmod(C3(state = 'A'), PARP(state='U'), PARP(bf = None, state='C'), kc3parpf, kc3parpr, kc3parpc)
twostepmod(C6(state='A'), C8(state='pro'), C8(bf = None, state = 'A'), kc6c8f, kc6c8r, kc6c8c)


#Initial states
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

# Fig 4B from Albeck observes these, normalizes and inverts them
Observe('Bid',   Bid(bf=None, state='U'))
Observe('PARP',  PARP(bf=None, state='U'))
Observe('Smac',  Smac(bf=None, state='mito'))
# This is what *should* be observed???
Observe('tBid',  Bid(state='T'))
Observe('cPARP', PARP(state='C'))
Observe('cSmac', Smac(state='cyto'))


# ####


# if __name__ == '__main__':
#     from pysb.bng import generate_network_code
#     from pysb.tools.export_bng import run as run_export
#     print run_export(model)
#     print ""
#     print "begin actions"
#     print "  generate_network({overwrite=>1});"
#     print "  simulate_ode({t_end=>21600,n_steps=>360});" # 6 hours, 1-minute steps
#     print "end actions"
