from pysb import *
from pysbhelperfuncs import *

# Instantiate model, default name is "model"
Model()

Parameter('ec_size', 1.0e6)   # 1.0e6 um^3 = 1 pL
Parameter('cytoM_size', 483.6 * .0030) # plasma SA (6.22um radius for a 1e3 um^3 cell) * membrane thickness ~3.0nm
Parameter('cyto_size', 1.0e3) # 1.0e3 um^3 --> size of HeLa. Range is 760-2730 um^3 (ref)
Parameter('mito_size', 70.0)  # mitochondria is ~7% of cytoplasm (ref)
Parameter('mitoM_size', 82.14 * .0042) # mito SA (2.55um radius) x 'brane thicknes ~4.2 nm JPC-B(2009)113-11p3413
Compartment('ec', dimension = 3, size = ec_size, parent = None)    # extra cellular compartment
Compartment('cyM', dimension = 2, size = cytoM_size, parent = ec) # cytoplasmic membrane
Compartment('cy', dimension = 3, size = cyto_size, parent = cyM)  # cytoplasm
Compartment('mitM', dimension = 2, size = mitoM_size, parent = cy) # mitochondrial membrane
Compartment('mit', dimension = 3, size = mito_size, parent = mitM) # mitochondrion

Monomer('Trail', ['b', 's1', 's2'])  # TRAIL monomer, start with pre-trimerized ligand
Monomer('DR', ['bl', 'bf', 's1', 's2', 'T'], {'T':['4','5']})    # Death receptor 4 or 5. bl: ligand binding site, bf: fadd binding site
Monomer('Fadd', ['bx', 'bc'])   # FADD bx: binding to complex, bc: binding to caspase 8
Monomer('flip', ['bf']) # flip
Monomer('C8', ['bc', 'bf', 'state'], {'state':['pro', 'A']}) # Csp 8, states: pro, active
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

# MONOMERS FROM EARM2
# Parameters and Modules 
# ===============================
from disc2_parms import parameter_dict as kd 
import earm_2_emb_modules # Must be called after the Monomers and Parameters are defined

# Trail binding to DR
# Another way to do this:
# bla= TTrail_unbound.copy()
# bla.monomer_patterns[0].site_conditions['bf']=1

# aliases for easier rule ennumeration
# -------------------------------------
# Trimerized Trail unbound/bound aliases
TTrail_U  = MatchOnce(Trail(b=None, s1=1, s2=2) % Trail(b=None, s1=2, s2=3) % Trail(b=None, s1=3, s2=1))
TTrail_B1 = MatchOnce(Trail(b=4,    s1=1, s2=2) % Trail(b=None, s1=2, s2=3) % Trail(b=None, s1=3, s2=1))
TTrail_B2 = MatchOnce(Trail(b=5,    s1=1, s2=2) % Trail(b=6,    s1=2, s2=3) % Trail(b=None, s1=3, s2=1)) 
TTrail_B3 = MatchOnce(Trail(b=4,    s1=1, s2=2) % Trail(b=5,    s1=2, s2=3) % Trail(b=6,    s1=3, s2=1))
TTrail_B3N = Trail(b=4,    s1=1, s2=2) % Trail(b=5,    s1=2, s2=3) % Trail(b=6,    s1=3, s2=1)

# DR monomer, dimer, trimer aliases:
# ----------------------------------
DR_mono_U = DR(bl=None, bf=None, s1=None, s2=None)
DR_dim_U  = DR(bl=None, bf=None, s1=1,    s2=None) % DR(bl=None, bf=None, s1=None, s2=1)
DR_trim_U = MatchOnce(DR(bl=None, bf=None, s1=1,    s2=2) % DR(bl=None, bf=None, s1=2,    s2=3) % DR(bl=None, bf=None, s1=3, s2=1))
DR_mono_B = DR(bl=4,    bf=None, s1=None, s2=None)
DR_dim_B  = MatchOnce(DR(bl=5,    bf=None, s1=4,    s2=None) % DR(bl=6, bf=None, s1=None, s2=4))
DR_trim_B = MatchOnce(DR(bl=4,    bf=None, s1=7,    s2=8   ) % DR(bl=5, bf=None, s1=8,    s2=9) % DR(bl=6,    bf=None, s1=9, s2=7))
DR_trim_BN =  DR(bl=4, bf=None, s1=7, s2=8  ) % DR(bl=5, bf=None, s1=8, s2=9) % DR(bl=6, bf=None, s1=9, s2=7)
DR_trim_BB1 = DR(bl=4,   bf=10, s1=7, s2=8  ) % DR(bl=5, bf=None, s1=8, s2=9) % DR(bl=6, bf=None,   s1=9, s2=7)
DR_trim_BB2 = DR(bl=4,   bf=10, s1=7, s2=8  ) % DR(bl=5, bf=11,   s1=8, s2=9) % DR(bl=6, bf=None, s1=9, s2=7)
DR_trim_BB3 = DR(bl=4,   bf=10, s1=7, s2=8  ) % DR(bl=5, bf=11,   s1=8, s2=9) % DR(bl=6, bf=12,   s1=9, s2=7)

# DR multimers
# ------------
ringp_assembly(DR(bl = None, bf=None, T='4'), 3, kd['DR4_RINGP'])
ringp_assembly(DR(bl = None, bf=None, T='5'), 3, kd['DR5_RINGP'])

# Trail binding to DR rules:
# --------------------------
# Paper Carlos R. Reis, Robbert H. Cool rates 2011
Rule('TRAIL_DRmono', TTrail_U + DR_mono_U <> TTrail_B1 % DR_mono_B, *kd['TT_DRmono'])
Rule('TRAIL_DRdim',  TTrail_U + DR_dim_U  <> TTrail_B2 % DR_dim_B,  *kd['TT_DRdim'])
Rule('TRAIL_DRtrim', TTrail_U + DR_trim_U <> TTrail_B3 % DR_trim_B, *kd['TT_DRtrim'])
    
# Ligand Death-receptor complex: 
# ------------------------------
LDRC =   TTrail_B3N % DR_trim_BN
LDRC_B1 = TTrail_B3N % DR_trim_BB1
LDRC_B2 = TTrail_B3N % DR_trim_BB2
LDRC_B3 = TTrail_B3N % DR_trim_BB3

# FIXME: Can these all be replaced with simple-bind calls?
# Fadd binding LDRC rule
# -------------------------------------------------------
Rule("Fadd_LDRC_cplx1", LDRC + Fadd(bx = None, bc=None) <> LDRC_B1 % Fadd(bx=10, bc=None), *kd['Fadd_LDRC'])
Rule("Fadd_LDRC_cplx2", LDRC_B1 % Fadd(bx=10, bc=None) + Fadd(bx = None, bc=None) <>
                        LDRC_B2 % Fadd(bx=10, bc=None) % Fadd(bx=11, bc=None), *kd['Fadd_LDRC'])
Rule("Fadd_LDRC_cplx3", LDRC_B2 % Fadd(bx=10, bc=None) % Fadd(bx=11, bc=None) + Fadd(bx = None, bc=None) <>
                        LDRC_B3 % Fadd(bx=10, bc=None) % Fadd(bx=11, bc=None) % Fadd(bx=12, bc=None), *kd['Fadd_LDRC'])

# Caspase 8 binding to Fadd, creates DISC
# ---------------------------------------
Rule("pC8_fadd_b", Fadd(bc=None, bx=ANY) + C8(bf=None, state='pro') <> Fadd(bc=1, bx=ANY) % C8(bf=1, state='pro'), *kd['pC8_fadd_b'])

# Caspase 8 activation within the same DISC
# -----------------------------------------------------------------
Rule("pC8_dim_act_s", C8(bc=None, state='pro') % C8(bc=None, state='pro') >> C8(bc=1, bf=None, state='A') % C8(bc=1, bf=None, state='A'), *kd['pC8_dim_act_s'])

# Caspase 8 activation across separate DISC
# -----------------------------------------
# FIXME: Keep in mind but do not include until JR gives a solid reference
# Rule("pC8_dim_act_o", C8(bf=ANY) + C8(bf=ANY) >> C8(bc=1, bf=None, state='A') % C8(bc=1, bf=None, state='A'), *kd['pC8_dim_act_o'])

# Truncation of Bid
# FIXME: use catalysis function here
# ----------------------------------
Rule("C8_bid_cplx", C8(bc=1, bf=None, state='A') % C8(bc=1, bf=None, state='A') + Bid(bf = None, state='U') <>
     C8(bc=1, bf=2, state='A') % C8(bc=1, bf=None, state='A') % Bid(bf = 2, state='U'), kd['C8_BID'][0], kd['C8_BID'][1])
Rule("C8_cplx_act", C8(bc=1, bf=2, state='A') % C8(bc=1, bf=None, state='A') % Bid(bf = 2, state='U') >>
     C8(bc=1, bf=None, state='A') % C8(bc=1, bf=None, state='A') + Bid(bf = None, state='T'), kd['C8_BID'][2])

# flip takes a C8 spot in DISC
# ----------------------------
# FIXME: THIS SHOULD NOT GENERATE Fadd:flip where fadd is not bound to a complex...?
Rule("flip_fadd_cplx", Fadd(bx=1, bc=None) + flip(bf=None) <>  Fadd(bx=1, bc=2) % flip(bf=2), *kd['DISC_FLIP'])

# BAR inhibits C8
# ---------------
#Rule("BAR_C8_cplx", C8(bc=1, bf=None, state='A') % C8(bc=1, bf=None, state='A') + BAR(bf = None) <>
#     C8(bc=1, bf=2, state='A') % C8(bc=1, bf=None, state='A') % BAR(bf = 2), kd['BAR_C8'][0], kd['BAR_C8'][1])

#
# Import necessary modules
# ========================
# Generate the Receptor to Bid section from the EARM 1.0 module
earm_2_emb_modules.bid_to_momp(model, kd)
# Generate the Pore to MOMP section from the EARM 1.0 module
earm_2_emb_modules.pore_to_parp(model, kd)


# ===================================
# 'ec'   # extra cellular compartment
# 'cyM'  # cytoplasmic membrane
# 'cy'   # cytoplasm
# 'mitM' # mitochondrial membrane
# 'mit'  # mitochondrion
# 
# Initial non-zero species
# ===================================
Initial((Trail(b=None, s1=1, s2=2) % Trail(b=None, s1=2, s2=3) % Trail(b=None, s1=3, s2=1)) ** ec, Trail_0)
Initial(DR(bl=None, bf=None, s1=None, s2=None, T='4') ** cyM, DR4_0)
Initial(DR(bl=None, bf=None, s1=None, s2=None, T='5') ** cyM, DR5_0)
Initial(Fadd(bx=None, bc=None) ** cyM, Fadd_0)
Initial(flip(bf=None) ** cyM, flip_0)
Initial(C8(bc=None, bf=None, state='pro') ** cy, C8_0)
Initial(BAR(bf=None) ** cy, BAR_0)
Initial(Bid(bf=None, state='U') ** cy, Bid_0)
Initial(Bax(bf=None, s1=None, s2=None, state='C') ** cy, Bax_0)
Initial(Bak(bf=None, s1=None, s2=None, state='M') ** mitM, Bak_0)
Initial(Bcl2(bf=None) ** mitM, Bcl2_0)
Initial(BclxL (bf=None, state='C') ** cy, BclxL_0)
Initial(Mcl1(bf=None) ** mitM, Mcl1_0)
Initial(Bad(bf=None, state='C') ** cy, Bad_0) 
Initial(NOXA(bf=None) ** mitM, NOXA_0)
Initial(CytoC(bf=None, state='M') ** mit, CytoC_0)
Initial(Smac(bf=None, state='M') ** mit, Smac_0)
Initial(Apaf(bf=None, state='I') ** cy, Apaf_0)
Initial(C3(bf=None, state='pro') ** cy, C3_0)
Initial(C6(bf=None, state='pro') ** cy, C6_0)
Initial(C9(bf=None) ** cy, C9_0)
Initial(PARP(bf=None, state='U') ** cy, PARP_0)
Initial(XIAP(bf=None) ** cy, XIAP_0)

# Observables
# ===========
# Fig 4B from Albeck observes these, normalizes and inverts them
# Observe('Bid',   Bid(bf=None, state='U'))
# Observe('PARP',  PARP(bf=None, state='U'))
# Observe('Smac',  Smac(bf=None, state='mito'))
# # This is what *should* be observed???
Observe('mBid',  Bid(state='M') ** mitM)
Observe('cSmac', Smac(state='A') ** cy)
Observe('cPARP', PARP(state='C') ** cy)
