from pysb import *
from pysbhelperfuncs import *

# Instantiate model, default name is "model"
Model()

#Parameter('ec_size', 1.0e6)   # 1.0e6 um^3 = 1 pL
#Parameter('cytoM_size', 483.6 * .0030) # plasma SA (6.22um radius for a 1e3 um^3 cell) * membrane thickness ~3.0nm
#Parameter('cyto_size', 1.0e3) # 1.0e3 um^3 --> size of HeLa. Range is 760-2730 um^3 (ref)
#Parameter('mito_size', 70.0)  # mitochondria is ~7% of citoplasm (ref)
#Parameter('mitoM_size', 82.14 * .0042) # mito SA (2.55um radius) x 'brane thicknes ~4.2 nm JPC-B(2009)113-11p3413
#Compartment('ec', dimension = 3, size = ec_size, parent = None)    # extra cellular compartment
#Compartment('cytM', dimension = 2, size = cytoM_size, parent = ec) # cytoplasmic membrane
#Compartment('cyt', dimension = 3, size = cyto_size, parent = cyM)  # cytoplasm
#Compartment('mitM', dimension = 2, size = mitoM_size, parent = cy) # mitochondrial membrane
#Compartment('mit', dimension = 3, size = mito_size, parent = mitM) # mitochondrion

#Monomer('C8', ['bf', 'state'], {'state':['pro', 'A']}) # Csp 8, states: pro, active
Monomer('Trail', ['b', 's1', 's2'])  # TRAIL monomer, start with pre-trimerized ligand
Monomer('DR', ['bl', 'bf', 's1', 's2', 'T'], {'T':['4','5']})    # Death receptor 4 or 5. bl: ligand binding site, bf: fadd binding site
Monomer('Fadd', ['bx', 'bc'])   # FADD bx: binding to complex, bc: binding to caspase 8
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

# DR monomer, dimer, trimer aliases:
# ----------------------------------
DR_mono_U = DR(bl=None, bf=None, s1=None, s2=None, T=ANY)
DR_dim_U  = DR(bl=None, bf=None, s1=1,    s2=None, T=ANY) % DR(bl=None, bf=None, s1=None, s2=1, T=ANY)
DR_trim_U = DR(bl=None, bf=None, s1=1,    s2=2,    T=ANY) % DR(bl=None, bf=None, s1=2,    s2=3, T=ANY) % DR(bl=None, bf=None, s1=3, s2=1, T=ANY)
DR_mono_B = DR(bl=4,    bf=None, s1=None, s2=None, T=ANY)
DR_dim_B  = DR(bl=5,    bf=None, s1=4,    s2=None, T=ANY) % DR(bl=6,    bf=None, s1=None, s2=4, T=ANY)
DR_trim_B = DR(bl=4,    bf=None, s1=1,    s2=2,    T=ANY) % DR(bl=5,    bf=None, s1=2,    s2=3, T=ANY) % DR(bl=6,    bf=None, s1=3, s2=1, T=ANY)

# Trail binding to DR rules:
# --------------------------
Rule('TRAIL_DRmono', TTrail_U + DR_mono_U <> TTrail_B1 % DR_mono_B, *kd['TT_DRmono'])
Rule('TRAIL_DRdim',  TTrail_U + DR_dim_U  <> TTrail_B2 % DR_dim_B,  *kd['TT_DRdim'])
Rule('TRAIL_DRtrim', TTrail_U + DR_trim_U <> TTrail_B3 % DR_trim_B, *kd['TT_DRtrim'])
    
# DR multimers
# ------------
ringp_assembly(DR(bf=None, T='4'), 3, kd['DR4_RINGP'])
ringp_assembly(DR(bf=None, T='5'), 3, kd['DR5_RINGP'])

# DR DISC aliases:
# ----------------
LDRC = DR_trim_B % TTrail_B3

LDRC_F = LDRC % Fadd(bx = 7)
LDRC_F.monomer_patterns[0].site_conditions['bf'] = 7

# FIXME: Can these all be replaced with simple-bind calls?
# Fadd binding LDRC rule
# ----------------------
# This should create a species which binds a Fadd to LDRC
Rule("Fadd_LDRC", LDRC + Fadd(bx = None) <> LDRC_F, kf, kr)

# Caspase 8 binding to Fadd, creates DISC
# ---------------------------------------
Rule("pC8_fadd_b", Fadd(bx=ANY, bc=None) + C8(bf=None, state='pro') <> Fadd(bx=ANY, bc=1) % C8(bf=1, state='pro'), kf, kr)

# Caspase 8 activation within the same DISC
# ------------------------------------------
Rule("pC8_dim_act_s", C8(bf=ANY) % C8(bf=ANY) >> C8(bc=1, bf=None, state='act') % C8(bc=1, bf=None, state='act'))

# Caspase 8 activation across separate DISC
# -----------------------------------------
Rule("pC8_dim_act_o", C8(bf=ANY) + C8(bf=ANY) >> C8(bc=1, bf=None, state='act') % C8(bc=1, bf=None, state='act'))

# Truncation of Bid
# ------------------
Rule("C8_bid_cplx", C8(bc=1, bf=None, state='act') % C8(bc=1, bf=None, state='act') + Bid(bf = None, state='U') <>
     C8(bc=1, bf=2, state='act') % C8(bc=1, bf=None, state='act') + Bid(bf = 2, state='U'))
Rule("C8_cplx_act", C8(bc=1, bf=2, state='act') % C8(bc=1, bf=None, state='act') + Bid(bf = 2, state='U') >>
     C8(bc=1, bf=2, state='act') % C8(bc=1, bf=None, state='act') + Bid(bf = None, state='T'))

# flip takes a C8 spot in DISC
# ----------------------------
Rule("flip_fadd_b", Fadd(bx=ANY, bc=None) + C8(bf=None, state='pro') <> Fadd(bx=ANY, bc=1) % C8(bf=1, state='pro'), *kd['pC8_fadd_b'])

simple_bind(BAR(), C8(state='A'), kd['BAR_C8'])


#
# Import necessary modules
# ========================
# Generate the Receptor to Bid section from the EARM 1.0 module
earm_2_emb_modules.bid_to_momp(model, kd)
# Generate the Pore to MOMP section from the EARM 1.0 module
earm_2_emb_modules.pore_to_parp(model, kd)


# Initial non-zero species
# ========================
Initial(Trail(bf=None, s1=1, s2=2) % Trail(bf=None, s1=2, s2=3) % Trail(bf=None, s1=3, s2=1), Trail_0)
Initial(R(bf=None), R_0)
