from pysb import *

"""
Model drawn from:

Chun Chen, Jun Cui, Haizhu Lu, Rui Wang, Shuai Zhang, and Pingping Shen (2007)
"Modeling of the role of a Bax-activation switch in the mitochondrial apoptosis
decision," Biophys J 92:4304-4315
"""
Model()

# Species

Monomer('tBid', ['bf'])
Monomer('Bax', ['bf', 'state'], {'state':['C', 'A']})
Monomer('Bcl2', ['bf'])
Monomer('Pore', ['bf'])

# Reactions 
one_step_mod(tBid(bf=None), Bax(bf=None, state='C'), Bax(bf=None, state='A', k1)
Rule('Bax_Inactivation', Bax(bf=None, state='A') >> Bax(bf=None, state='C'), k2)
simple_bind(Bax(state='A'), Bcl2(), [k3, k4]) # Replace these with bind table?
simple_bind(tBid(), Bcl2(), [k5, k6])
Rule('Bcl2_Exchange?',
     Bax(bf=None, state='A') + tBid(b=1) % Bcl2(b=1) <>
     tBid(b=None) + Bax(bf=1, state='C') % Bcl2(b=1), k7, k8)
#Rule('Bax_Tetramerization', AcBax(b=None) + AcBax(b=None) + AcBax(b=None) + AcBax(b=None) <> Bax4(b=None), k9, k10)

# Observables
Observe('Act',  Act())
#Observe('Act_free', Act(b=None))
Observe('InBax_tot', Bax())
#Observe('Bax_free', Bax(b=None))
#Observe('AcBax_free',  AcBax(b=None))
Observe('AcBax_tot', Bax(s))
#Observe('AcBax_Bcl2', AcBax(b=1)%Bcl2(b=1))
Observe('Bcl2_free', Bcl2(b=None))
#Observe('Bax4', Bax4())


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
"""
# Starting concentrations at the mito
Act_mito = 0       # Activator 100% cytosolic
Bax_mito = 0.2
Bcl2_mito = 0.1
# Starting concentrations in the cytosol
Act_cyto = 0.1
Bax_cyto = 0.2
Bcl2_cyto = 0      # Bcl2 100% mitochondrial
# Constants to scale the amount of translocation from cyto to mito
F1_Bax =  0      # Fraction of Bax relocalized
F2_Act =  0.01      # Fraction of Activator relocalized
F3_Bcl2 = 0      # Fraction of Bcl2 relocalized

# Now, set the initial conditions at the mito
Parameter('Bax_0', Bax_mito + (Bax_cyto * F1_Bax))
Parameter('Act_0', Act_mito + (Act_cyto * F2_Act))
Parameter('Bcl2_0', Bcl2_mito + (Bcl2_cyto * F3_Bcl2))
# Rate parameters
Parameter('k1', 0.5)   # uM-1 s-1
Parameter('k2', 0.1)   # s-1
Parameter('k3', 2)     # um-1 s-1
Parameter('k4', 0.001)
Parameter('k5', 3)     # um-1 s-1
Parameter('k6', 0.04)  # s-1
Parameter('k7', 2)     # um-1 s-1
Parameter('k8', 0)
Parameter('k9', 2)     # um-1 s-1
Parameter('k10', 0)
"""

