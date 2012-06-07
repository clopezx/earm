from pysb import *

# Albeck JG, Burke JM, Spencer SL, Lauffenburger DA, Sorger PK, 2008
# Modeling a Snap-Action, Variable-Delay Switch Controlling Extrinsic
# Cell Death. PLoS Biol 6(12): e299. doi:10.1371/journal.pbio.0060299
#

Model()

def catalyze(enz, sub, prod, kf, kr, kc):
    """2-step catalytic process"""
    r1_name = 'bind_%s_%s' % (sub.name, enz.name)
    r2_name = 'produce_%s_via_%s' % (prod.name, enz.name)
    E = enz(b=None)
    S = sub(b=None)
    ES = enz(b=1) % sub(b=1)
    P = prod(b=None)
    Rule(r1_name, E + S <> ES, kf, kr)
    Rule(r2_name, ES >> E + P, kc)

def catalyze_convert(s1, s2, p, kf, kr, kc):
    """2-step catalytic-type process, but the "catalyst" is effectively consumed"""
    r1_name = 'bind_%s_%s' % (s1.name, s2.name)
    r2_name = 'produce_%s' % p.name
    A = s1(b=None)
    B = s2(b=None)
    AB = s1(b=1) % s2(b=1)
    C = p(b=None)
    Rule(r1_name, A + B <> AB, kf, kr)
    Rule(r2_name, AB >> C, kc)


def inhibit(targ, inh, kf, kr):
    """inhibition by complexation/sequestration"""
    r_name = 'inhibit_%s_by_%s' % (targ.name, inh.name)
    T = targ(b=None)
    I = inh(b=None)
    TI = targ(b=1) % inh(b=1)
    Rule(r_name, T + I <> TI, kf, kr)


# Species
Monomer('Act', ['b'])
Monomer('Bax', ['b'])
Monomer('AcBax', ['b'])
Monomer('Bcl2', ['b'])
Monomer('Bax4', ['b'])

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

# Reaction rules
Rule('Bax_Activation', Act(b=None) + Bax(b=None) >> Act(b=None) + AcBax(b=None), k1)
Rule('Bax_Inactivation', AcBax(b=None) >> Bax(b=None), k2)
Rule('Bcl2_Inhibits_Bax', AcBax(b=None) + Bcl2(b=None) <> AcBax(b=1) % Bcl2(b=1), k3, k4)
Rule('Bcl2_Inhibits_Activator', Act(b=None) + Bcl2(b=None) <> Act(b=1) % Bcl2(b=1), k5, k6)
Rule('Bcl2_Exchange?', AcBax(b=None) + Act(b=1)%Bcl2(b=1) <> Act(b=None) + AcBax(b=1)%Bcl2(b=1), k7, k8)
#Rule('Bax_Tetramerization', AcBax(b=None) + AcBax(b=None) + AcBax(b=None) + AcBax(b=None) <> Bax4(b=None), k9, k10)

# Observables
Observe('Act',  Act())
#Observe('Act_free', Act(b=None))
Observe('InBax_tot', Bax())
#Observe('Bax_free', Bax(b=None))
#Observe('AcBax_free',  AcBax(b=None))
Observe('AcBax_tot', AcBax())
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

if __name__ == '__main__':
    from pysb.generator.bng import BngGenerator
    gen = BngGenerator(model)
    print gen.get_content()
    print ""
    print "begin actions"
    print "  generate_network({overwrite=>1});"
    print "  simulate_ode({t_end=>21600,n_steps=>360});" # 6 hours, 1-minute steps
    print "end actions"
