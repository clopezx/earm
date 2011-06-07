from pysb import *
from pysbhelperfuncs import *

L = Monomer('L', ['bf']) # Ligand
R = Monomer('R', ['bf']) # Receptor
DISC = Monomer('DISC', ['bf']) # DISC
flip = Monomer('flip', ['bf']) # flip
C8 = Monomer('C8', ['bf', 'state'], {'state':['pro', 'A']}) # Csp 8, states: pro, active
BAR = Monomer('BAR', ['bf']) # BAR

# RECEPTOR TO tBID
# =====================
# tBID Activation Rules
# ---------------------
#        L + R <--> L:R --> DISC
#        pC8 + DISC <--> DISC:pC8 --> C8 + DISC
#        Bid + C8 <--> Bid:C8 --> tBid + C8
# ---------------------
twostepconv(L(), R(), DISC(), klrf, klrr, klrc)
twostepmod(DISC(), C8(state='pro'), C8(bf = None, state='A'), kdiscc8f, kdiscc8r, kdiscc8c)
twostepmod(C8(state='A'), Bid(state='U'), Bid(state='T'), kc8bidf, kc8bidr, kc8bidc)
# ---------------------
# Inhibition Rules
# ---------------------
#        flip + DISC <-->  flip:DISC  
#        C8 + BAR <--> BAR:C8 CSPS
# ---------------------
simplebind(DISC(), flip(), kflipdiscf, kflipdiscr)
simplebind(BAR(), C8(state='A'), kbarc8f, kbarc8r)

