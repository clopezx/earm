from pysb import *
from pysbhelperfuncs import *

Parameter('ec_size', 1.0e6)   # 1.0e6 um^3 = 1 pL
Parameter('cyto_size', 1.0e3) # 1.0e3 um^3 --> size of HeLa. Range is 760-2730 um^3 (ref)
Parameter('mito_size', 70.0)  # mitochondria is ~7% of citoplasm (ref)
Parameter('mitoM_size', 82.14 * .0042) # mitochondria SA (2.55um radius) * membrane thicknes ~4.2 nm
                                       # J. Phys. Chem. B, 2009, 113 (11), pp 3413–3422 DOI: 10.1021/jp8077369
