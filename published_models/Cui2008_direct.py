from pysb import *
from local_macros import *
from shen_direct_module import direct_module

Model()

declare_monomers()

# MODEL TOPOLOGY
direct_module()

synthesize_degrade_table(
   [[Bax(state='C'),                   0.06,  0.001],
    [Bax(state='A'),                   None,  0.001],                     
    [tBid,                            0.001,  0.001],
    [Bcl2,                             0.03,  0.001],
    [tBid(b=1) % Bcl2(b=1),            None,  0.005],
    [Bax(b=1, state='C') % Bcl2(b=1),  None,  0.005],
    [Bad,                             0.001,  0.001],
    [Bad(b=1) % Bcl2(b=1),             None,  0.005],
    [Pore,                             None, 0.0005]])
