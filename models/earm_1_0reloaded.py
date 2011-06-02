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
import earm_1_0_parms

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
Monomer('Bid', ['bf'], 'state'], {'state':['U', 'T']}) # Bid, states: untruncated, truncated
Monomer('Bcl2', ['bf', 'state'], {'state':['cyto', 'mito'}) # Bcl2, states: cytoplasm, mitochondria
Monomer('Bax', ['bf', 'state'], {'state':['I', 'A', 'M']}) # Bax, states: inactive, active, mitochondrial
Monomer('Bax2', ['bf'])
Monomer('Bax4', ['bf'])
Monomer('MitoP', ['bf', 'state'],{'state':['U', 'A']})
Monomer('CytoC', ['bf', 'state'], {'state':['mito', 'A', 'cyto']}
Monomer('Smac', ['bf', 'state'], {'state':['mito', 'A', 'cyto']}
Monomer('C9', ['bf'])
Monomer('Apaf', ['bf', 'state'], {'state':['I', 'A']}
Monomer('Apop', ['bf'])
Monomer('XIAP', ['bf'])

# Rules
# Ligand + Receptor <--> L:R --> DISC
twostepconv(L(), R(), DISC(), klrf, klrr, klrc)

# flip + DISC <-->  flip:DISC  
simplebind(DISC(), flip(), kflipdiscf, kflipdiscr)

# pC8 + DISC <--> DISC:pC8 --> C8 + DISC
twostepact(DISC(), C8(state='pro'), C8(bf = None, state='A'),
           kdiscc8f, kdiscc8r, kdiscc8c)

# C8 + BAR <--> BAR:C8 
simplebind(BAR(), C8(state='A'), kbarc8f, kbarc8r)

# pC3 + C8 <--> pC3:C8 --> C3 + C8
twostepact(C8(state='A'), C3(state='pro'), C3(bf = None, state='A'),
           kc8c3f, kc8c3r, kc8c3c)

# pC6 + C3 <--> pC6:C3 --> C6 + C3
twostepact(C3(state='A'), C6(state='pro'), C6(bf = None, state='A'),
           kc3c6f, kc3c6r, kc3c6c)

# pC8 + C6 <--> pC8:C6 --> C8 + C6
twostepact(C6(state='A'), C8(state='pro'), C8(bf = None, state = 'A'),
           kc6c8f, kc6c8r, kc8c6c)

# XIAP + C3 <--> XIAP:C3 --> XIAP + C3_U
twostepact(XIAP, C3(state = 'A'), C3(bf = None, state = 'ub'),
           kxiapc3f, kxiapc3r, kxiapc3c)

# PARP + C3 <--> PARP:C3 --> CPARP + C3
twostepact(C3(state = 'A'), PARP(state='U'), PARP(bf = None, state='C')

# Bid + C8 <--> Bid:C8 --> tBid + C8
twostepact(C8(state='A'), Bid(state='U'), Bid(state='T'),

# tBid + Bcl2c <-->  tBid:Bcl2c INHIBITION IN CYTO ONLY
simplebind(Bid(state='T'), Bcl2(state='cyto'), kbidbcl2f, kbidbcl2r)

# Bax + tBid <--> Bax:tBid --> aBax + tBid 
twostateact(Bid(state = 'T'), Bax(state='I'), Bax(bf = None, state = 'A'),
            kbidbaxf, kbidbaxr, kbidbaxc)

# aBax <-->  MBax 
Rule('baxCtoM', Bax(bf = None, state = 'A') <> Bax(bf=None, state = 'M'),
     kbaxcbaxmf, kbaxcbaxmr)

# MBax + Bcl2 <-->  MBax:Bcl2  
simplebind(Bax(state='M'), Bcl2(state='mito'), kbaxMbcl2Mf, kbaxMbcl2Mr)

# MBax + MBax <-->  Bax2
simpledimer(Bax(state='M'), Bax2(bf=None), kbaxdimf, kbaxdimr)
###Rule('dimerize_MBax_to_Bax2', MBax(b=None) + MBax(b=None) <> Bax2(b=None), kf15, kr15)

# Bax2 + Bcl2 <-->  Bax2:Bcl2  
simplebind(Bax2(), Bcl2(state='mito'), kbax2Mbcl2Mf, kbax2Mbcl2Mr)

# Bax2 + Bax2 <-->  Bax4
simpledimer(Bax2(), Bax4(bf = None), kbaxtetf, kbaxtetr)
## Rule('dimerize_Bax2_to_Bax4', Bax2(b=None) + Bax2(b=None) <> Bax4(b=None), kf17, kr17)

# Bax4 + Bcl2 <-->  Bax4:Bcl2  
simplebind(Bax4(), Bcl2(state='mito'), kbax4Mbcl2Mf, kbax4Mbcl2Mr)

# Bax4 + Mito <-->  Bax4:Mito -->  AMito  
twostepconv(Bax4(), MitoP(state='U'), MitoP(bf = None, state='A'),
            kbax4poref,kbax4porer,kbax4porec) 

# AMito + mCytoC <-->  AMito:mCytoC --> AMito + ACytoC  
twostepact(MitoP(state='A'), CytoC(state='mito'), CytoC(bf = None, state='A'),
           kmitopcytocMf, kmitopcytocMr, kmitopcytocMc)

# AMito + mSmac <-->  AMito:mSmac --> AMito + ASmac  
twostepact(MitoP(state='A'), Smac(state='mito'), Smac(bf = None, state='A'),
           kmitopsmacMf, kmitopsmacMr, kmitopsmacMc)

# ACytoC <-->  cCytoC
Rule('cytocCtoM', CytoC(bf = None, state = 'A') <> CytoC(bf=None, state = 'cyto'),
     kcytocMcytocCf, kcytocMcytocCr)

# Apaf + cCytoC <-->  Apaf:cCytoC --> aApaf + cCytoC
twostepact(CytoC(state='cyto'), Apaf(state='I'), Apaf(bf = None, state = 'A'),
           kcytocCapaff, kcytocCapafr, kcytocCapafc)

# aApaf + pC9 <-->  Apop
simpleconv(Apaf(state='A', C9(), Apop(bf=None),
                kapafc9f, kapafc9r)

# Apop + pC3 <-->  Apop:pC3 --> Apop + C3
twostepact(Apop(), C3(state='pro'), C3(bf = None, state='A'),
           kapopc3f, kapopc3r, kapopc3c)

# ASmac <-->  cSmac
Rule('SmacMtoSmacC', Smac(state='A') <> Smac(state='cyto'), ksmacMsmacCf, ksmacMsmacCr)

# Apop + XIAP <-->  Apop:XIAP  
simplebind(Apop(), XIAP(), kapopxiapf, kapopxiapr)

# cSmac + XIAP <-->  cSmac:XIAP  
simplebind(Smac(state='cyto'), XIAP(), ksmacxiapf, ksmacxiapr)


# Fig 4B
Observe('Bid',   Bid(b=None))
Observe('PARP',  PARP(b=None))
Observe('mSmac', mSmac(b=None))
# this is what I originally thought 4B was actually plotting
Observe('tBid',  tBid())
Observe('CPARP', CPARP())
Observe('cSmac', cSmac())



# generate initial conditions from _0 parameter naming convention
for m in model.monomers:
    ic_param = model.parameter('%s_0' % m.name)
    if ic_param is not None:
        sites = {}
        for s in m.sites:
            if s in m.site_states:
                sites[s] = m.site_states[s][0]
            else:
                sites[s] = None
        Initial(m(sites), ic_param)


####


if __name__ == '__main__':
    from pysb.bng import generate_network_code
    from pysb.tools.export_bng import run as run_export
    print run_export(model)
