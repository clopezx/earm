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
Monomer('pC9', ['bf'])
Monomer('Apaf', ['bf', 'state'], {'state':['I', 'A']}
Monomer('Apop', ['bf'])




# Ligand + Receptor <--> L:R --> DISC
twostepconv(L(), R(), DISC(), klrf, klrr, klrc)

# flip + DISC <-->  flip:DISC  
simplebind(DISC(), flip(), kflipdiscf, kflipdiscr)

# pC8 + DISC <--> DISC:pC8 --> C8 + DISC
twostepact(DISC(), C8(state='pro'), C8(state='A'),
           kdiscc8f, kdiscc8r, kdiscc8c)

# C8 + BAR <--> BAR:C8 
Parameter('kf4', 1e-06)
Parameter('kr4', 1e-03)
inhibit(BAR, C8, kf4, kr4)

# pC3 + C8 <--> pC3:C8 --> C3 + C8
Parameter('kf5', 1e-07)
Parameter('kr5', 1e-03)
Parameter('kc5', 1e+00)
catalyze(C8, pC3, C3, kf5, kr5, kc5)

# pC6 + C3 <--> pC6:C3 --> C6 + C3
Parameter('kf6', 1e-06)
Parameter('kr6', 1e-03)
Parameter('kc6', 1e+00)
catalyze(C3, pC6, C6, kf6, kr6, kc6)

# pC8 + C6 <--> pC8:C6 --> C8 + C6
Parameter('kf7', 3e-08)
Parameter('kr7', 1e-03)
Parameter('kc7', 1e+00)
catalyze(C6, pC8, C8, kf7, kr7, kc7)

# XIAP + C3 <--> XIAP:C3 --> XIAP + C3_U
Parameter('kf8', 2e-06)
Parameter('kr8', 1e-03)
Parameter('kc8', 1e-01)
catalyze(XIAP, C3, C3_U, kf8, kr8, kc8)

# PARP + C3 <--> PARP:C3 --> CPARP + C3
Parameter('kf9', 1e-06)
Parameter('kr9', 1e-02)
Parameter('kc9', 1e+00)
catalyze(C3, PARP, CPARP, kf9, kr9, kc9)

# Bid + C8 <--> Bid:C8 --> tBid + C8
Parameter('kf10', 1e-07)
Parameter('kr10', 1e-03)
Parameter('kc10', 1e+00)
catalyze(C8, Bid, tBid, kf10, kr10, kc10)

# tBid + Bcl2c <-->  tBid:Bcl2c  
Parameter('kf11', 1e-06)
Parameter('kr11', 1e-03)
inhibit(tBid, Bcl2c, kf11, kr11)

# Bax + tBid <--> Bax:tBid --> aBax + tBid 
Parameter('kf12', 1e-07)
Parameter('kr12', 1e-03)
Parameter('kc12', 1e+00)
catalyze(tBid, Bax, aBax, kf12, kr12, kc12)

# aBax <-->  MBax 
Parameter('kf13', transloc)
Parameter('kr13', transloc)
Rule('transloc_MBax_aBax', aBax(b=None) <> MBax(b=None), kf13, kr13)

# MBax + Bcl2 <-->  MBax:Bcl2  
Parameter('kf14', 1e-06/v)
Parameter('kr14', 1e-03)
inhibit(MBax, Bcl2, kf14, kr14)

# MBax + MBax <-->  Bax2
Parameter('kf15', 1e-06/v*2)
Parameter('kr15', 1e-03)
Rule('dimerize_MBax_to_Bax2', MBax(b=None) + MBax(b=None) <> Bax2(b=None), kf15, kr15)

# Bax2 + Bcl2 <-->  Bax2:Bcl2  
Parameter('kf16', 1e-06/v)
Parameter('kr16', 1e-03)
inhibit(Bax2, Bcl2, kf16, kr16)

# Bax2 + Bax2 <-->  Bax4
Parameter('kf17', 1e-06/v*2)
Parameter('kr17', 1e-03)
Rule('dimerize_Bax2_to_Bax4', Bax2(b=None) + Bax2(b=None) <> Bax4(b=None), kf17, kr17)

# Bax4 + Bcl2 <-->  Bax4:Bcl2  
Parameter('kf18', 1e-06/v)
Parameter('kr18', 1e-03)
inhibit(Bax4, Bcl2, kf18, kr18)

# Bax4 + Mito <-->  Bax4:Mito -->  AMito  
Parameter('kf19', 1e-06/v)
Parameter('kr19', 1e-03)
Parameter('kc19', 1e+00)
catalyze_convert(Bax4, Mito, AMito, kf19, kr19, kc19)

# AMito + mCytoC <-->  AMito:mCytoC --> AMito + ACytoC  
Parameter('kf20', 2e-06/v)
Parameter('kr20', 1e-03)
Parameter('kc20', 1e+01)
catalyze(AMito, mCytoC, ACytoC, kf20, kr20, kc20)

# AMito + mSmac <-->  AMito:mSmac --> AMito + ASmac  
Parameter('kf21', 2e-06/v)
Parameter('kr21', 1e-03)
Parameter('kc21', 1e+01)
catalyze(AMito, mSmac, ASmac, kf21, kr21, kc21)

# ACytoC <-->  cCytoC
Parameter('kf22', transloc)
Parameter('kr22', transloc)
Rule('transloc_cCytoC_ACytoC', ACytoC(b=None) <> cCytoC(b=None), kf22, kr22)

# Apaf + cCytoC <-->  Apaf:cCytoC --> aApaf + cCytoC
Parameter('kf23', 5e-07)
Parameter('kr23', 1e-03)
Parameter('kc23', 1e+00)
catalyze(cCytoC, Apaf, aApaf, kf23, kr23, kc23)

# aApaf + pC9 <-->  Apop
Parameter('kf24', 5e-08)
Parameter('kr24', 1e-03)
Rule('bind_aApaf_pC9_as_Apop', aApaf(b=None) + pC9(b=None) <> Apop(b=None), kf24, kr24)

# Apop + pC3 <-->  Apop:pC3 --> Apop + C3
Parameter('kf25', 5e-09)
Parameter('kr25', 1e-03)
Parameter('kc25', 1e+00)
catalyze(Apop, pC3, C3, kf25, kr25, kc25)

# ASmac <-->  cSmac
Parameter('kf26', transloc)
Parameter('kr26', transloc)
Rule('transloc_cSmac_ASmac', ASmac(b=None) <> cSmac(b=None), kf26, kr26)

# Apop + XIAP <-->  Apop:XIAP  
Parameter('kf27', 2e-06)
Parameter('kr27', 1e-03)
inhibit(Apop, XIAP, kf27, kr27)

# cSmac + XIAP <-->  cSmac:XIAP  
Parameter('kf28', 7e-06)
Parameter('kr28', 1e-03)
inhibit(cSmac, XIAP, kf28, kr28)



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
