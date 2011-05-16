import string

def catalyzeori(enz, sub, prod, kf, kr, kc):
    """2-step catalytic process"""
    r1_name = 'bind_%s_%s' % (sub.name, enz.name)
    r2_name = 'produce_%s_via_%s' % (prod.name, enz.name)
    E = enz(b=None)
    S = sub(b=None)
    ES = enz(b=1) % sub(b=1)
    P = prod(b=None)
    Rule(r1_name, E + S <> ES, kf, kr)
    Rule(r2_name, ES >> E + P, kc)

def catalyze(R0, R1, P0, kf, kr, kc):
    """Automation of the R0 + R1 <> R0:R1 >> R0+P0 two-step catalytic reaction
    this function assumes that there is a site named 'bc' (bind catalysis)
    hich it uses by default."""
    r1_name = 'bind_%s_%s' % (sub.name, enz.name)
    r2_name = 'produce_%s_via_%s' % (prod.name, enz.name)
    try:
        R0.site_conditions['bc'] is None
        R1.site_conditions['bc'] is None
    except KeyError:
        raise KeyError('Required site 'bc' not present in %s or %s as required')
    #build the intermediate from R0 and R1 using 'bc' as the binding site
    templist=[]
    for x,y in R0.site_conditions.iteritems():
        #test for state sites, these contain strings in the value
        if type(y) is str:
            templist.append("%s='%s' "%(x,y))
        else:
            templist.append("%s=%s "%(x,y))
    
    
    R0str=""
    for x,y in R1.site_conditions.iteritems():
        if type(y) is str:
            R1str += "%s='%s' "%(x,y)
        else:
            R1str += "%s=%s "%(x,y)
    #build the R0:R1 string and exec it
    
    




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


