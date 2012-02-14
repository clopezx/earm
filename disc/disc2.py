from pysb import *
from pysbhelperfuncs import *

Parameter('ec_size', 1.0e6)   # 1.0e6 um^3 = 1 pL
Parameter('cytoM_size', 483.6 * .0030) # plasma SA (6.22um radius for a 1e3 um^3 cell) * membrane thickness ~3.0nm
Parameter('cyto_size', 1.0e3) # 1.0e3 um^3 --> size of HeLa. Range is 760-2730 um^3 (ref)
Parameter('mito_size', 70.0)  # mitochondria is ~7% of citoplasm (ref)
Parameter('mitoM_size', 82.14 * .0042) # mitochondria SA (2.55um radius) * membrane thicknes ~4.2 nm
                                       # J. Phys. Chem. B, 2009, 113 (11), pp 3413–3422 DOI: 10.1021/jp8077369

Compartment('ec', dimension = 3, size = ec_size, parent = None)    # extra cellular compartment
Compartment('cytM', dimension = 2, size = cytoM_size, parent = ec) # cytoplasmic membrane
Compartment('cyt', dimension = 3, size = cyto_size, parent = cyM)  # cytoplasm
Compartment('mitM', dimension = 2, size = mitoM_size, parent = cy) # mitochondrial membrane
Compartment('mit', dimension = 3, size = mito_size, parent = mitM) # mitochondrion

Monomer('TRAIL', [b])  # TRAIL trimerized ligand (??? check)
Monomer('DR', [b], {'T':[4,5]})    # Death receptor 4
Monomer('Fadd', [b])   # FADD
Monomer('LDRC', [b])   # Ligand :: Death receptor complex
Monomer('LFDRC', [b])  # Ligand :: Fadd :: Death receptor
# MONOMERS FROM EARM2

# Parameters and Modules 
# ===============================
from earm_2_disc_parms import parameter_dict as kd 
import earm_2_disc_modules # Must be called after the Monomers and Parameters are defined

# Trail binding to DR
Rule('TRAIL_bind', TRAIL(b=None) + DR(b=None, T=any) <> TRAIL(b=1) % DR(b=1, T=any),
     kd['TRAIL_bind'])

*** CAN TRAIL BE BOUND INDIVIDUALLY TO EACH DR? OR TO JUST ONE?

# Trail assembly (use oligomerization function)
complex_assembly(DR(b=None, T=Any), 3, kd['DR_oligomer']) #*** CHECK REVERSIBLE

# Fadd assembly
Rule('Fadd_bind', Fadd(b=None) + [])

# C8 binding

# C8 activation within same DISC

# C8 cross activation across two DISC

# Bid truncated by C8

# C3 Activated by C8

# IC-RP Activated by C8




