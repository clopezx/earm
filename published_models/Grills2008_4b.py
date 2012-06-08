from pysb import *
from local_macros import *

Model()

declare_monomers()

two_state_equilibrium(Bax(state='C'), Bax(state='A'), [1, 1])

bind_table([[                   Bcl2,    Mcl1 ],
            [ Bax(state='C'),  (1,1),   (1,1) ]])


