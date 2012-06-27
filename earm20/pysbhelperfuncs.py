import inspect
from pysb import *
from pysb.core import SelfExporter, MonomerPattern, ComplexPattern

def del_rule(model, rname):
    """delete rules by name 'rname'"""
    idx = [i for i,r in enumerate(model.rules) if r.name == rname][0]
    model.rules.pop(idx)

def report_extra_parms(model):
    """report the parameters that are not involved in any rules"""
    pass


def alias_model_components(model=None):
    """
    """
    if model is None:
        model = SelfExporter.default_model
    caller_globals = inspect.currentframe().f_back.f_globals
    caller_globals.update(model.all_components())

def two_step_mod(Enz, Sub, Prod, klist,  site='bf'):
    """Automation of the Enz + Sub <> Enz:Sub >> Enz + Prod two-step catalytic reaction.
    This function assumes that there is a site named 'bf' (bind site for fxn)
    which it uses by default. This also assume Enz returns to its original state.
    In an aim for simplicity, site 'bf' need not be passed when calling the function by
    the reactants"""
    
    kf, kr, kc = klist #forward, reverse, catalytic
    
    # FIXME: this will fail if the argument passed is a Complex object. 

    r1_name = 'cplx_%s%s_%s%s' % (Sub.monomer.name, ''.join(filter(lambda a: a != None, Sub.site_conditions.values())),
                                     Enz.monomer.name, ''.join(filter(lambda a: a != None, Enz.site_conditions.values())))
    r2_name = 'diss_%s%s_via_%s%s' % (Prod.monomer.name, ''.join(filter(lambda a: a != None, Prod.site_conditions.values())),
                                         Enz.monomer.name, ''.join(filter(lambda a: a != None, Enz.site_conditions.values())))
    
    assert site in Enz.monomer.sites_dict, \
        "Required site %s not present in %s as required"%(site, Enz.monomer.name)
    assert site in Sub.monomer.sites_dict, \
        "Required site %s not present in %s as required"%(site, Sub.monomer.name)

    # make the intermediate complex components
    etmpdict = Enz.site_conditions.copy()
    stmpdict = Sub.site_conditions.copy()
    
    etmpdict[site] = 1
    stmpdict[site] = 1

    EnzCplx = Enz.monomer(etmpdict)
    SubCplx = Sub.monomer(stmpdict)

    # add the 'bf' site to the patterns
    Enz.site_conditions[site] = None
    Sub.site_conditions[site] = None

    # now that we have the complex elements formed we can write the first step rule
    Rule(r1_name, Enz + Sub <> EnzCplx % SubCplx, kf, kr)
    
    # and finally the rule for the catalytic transformation
    Rule(r2_name, EnzCplx % SubCplx >> Enz + Prod, kc)

def two_step_conv(Sub1, Sub2, Prod, klist, site='bf'):
    """Automation of the Sub1 + Sub2 <> Sub1:Sub2 >> Prod two-step reaction (i.e. dimerization).
    This function assumes that there is a site named 'bf' (bind site for fxn)
    which it uses by default. Site 'bf' need not be passed when calling the function."""

    kf, kr, kc = klist
    
    r1_name = 'cplx_%s_%s' % (Sub2.monomer.name, Sub1.monomer.name)

    #FIXME: this is a bit dirty but it fixes the problem when prod is a pattern
    if isinstance(Prod, MonomerPattern):
        r2_name = 'cplx_%s_via_%s__%s' % (Prod.monomer.name, Sub1.monomer.name, Sub2.monomer.name)
    elif isinstance(Prod, ComplexPattern):
        r2_name = 'cplx_%s_via_%s__%s' % (("_".join([sub.monomer.name for sub in Prod.monomer_patterns])),
                                          Sub1.monomer.name, Sub2.monomer.name)
    
    assert site in Sub1.monomer.sites_dict, \
        "Required site %s not present in %s as required"%(site, Sub1.monomer.name)
    assert site in Sub2.monomer.sites_dict, \
        "Required site %s not present in %s as required"%(site, Sub2.monomer.name)

    # make the intermediate complex components
    s1tmpdict = Sub1.site_conditions.copy()
    s2tmpdict = Sub2.site_conditions.copy()
    
    s1tmpdict[site] = 1
    s2tmpdict[site] = 1

    Sub1Cplx = Sub1.monomer(s1tmpdict)
    Sub2Cplx = Sub2.monomer(s2tmpdict)

    # add the site to the patterns
    Sub1.site_conditions[site] = None
    Sub2.site_conditions[site] = None

    # now that we have the complex elements formed we can write the first step rule
    Rule(r1_name, Sub1 + Sub2 <> Sub1Cplx % Sub2Cplx, kf, kr)
    
    # and finally the rule for the catalytic transformation
    Rule(r2_name, Sub1Cplx % Sub2Cplx >> Prod, kc)

def simple_dim(Sub, Prod, klist, site='bf'):
    """ Convert two Sub species into one Prod species:
    Sub + Sub <> Prod
    """
    kf, kr = klist
    r1_name = 'dimer_%s_to_%s'%(Sub.monomer.name, Prod.monomer.name)
    assert site in Sub.monomer.sites_dict, \
        "Required site %s not present in %s as required"%(site, Sub.monomer.name)

    # create the sites for the monomers
    Sub.site_conditions[site] = None

    # combine the monomers into a product step rule
    Rule(r1_name, Sub + Sub <> Prod, kf, kr)

def ringp_species(Subunit, size, site1 = 's1', site2 = 's2'):
    """
    Generate a single species representing a homomeric ringp, composed
    of <size> copies of <Subunit> bound together in a ring, with bonds
    formed between s1 of one unit and s2 of the next.
    """

    if size == 0:
        raise ValueError("size must be an integer greater than 0")
    if size == 1:
        Ringp = Subunit({site1: None, site2: None})
    elif size == 2:
        Ringp = Subunit({site1: 1, site2: None}) % Subunit({site1: None, site2: 1})
    else:
        Ringp = ComplexPattern([], None, match_once=True)
        for i in range(1, size+1):
            Ringp %= Subunit({site1: i, site2: i%size+1})
    return Ringp

def ringp_assembly(Subunit, size, rates):
    """
    Generate rules to chain identical MonomerPatterns <Subunit> into
    increasingly larger ringps of up to <size> units, using sites
    bh3 and d2 to bind the units to each other.
    """
    rules = []
    for i in range(2, size + 1):
        M = ringp_species(Subunit, 1)
        S1 = ringp_species(Subunit, i-1)
        S2 = ringp_species(Subunit, i)
        # ensure the rule name is different for monomers with differing states
        rname = Subunit.monomer.name + "_" + "_".join(["".join(map(str, j)) for j in Subunit.site_conditions.iteritems()]) + "_" + str(i) + "mer"
        #print rname
        rules.append(Rule(rname, M + S1 <> S2, *rates[i-2]))
    return rules

def ringp_transport(Subunit, Source, Dest, min_size, max_size, rates, site='bf'):
    """
    Generate rules to transport MonomerPattern <Source> to <Dest>
    through any of a series of ringps of at least <min_size> and at
    most <max_size> subunits, as defined by ringp_assembly. Uses site
    <site> on both Subunit and CargoSource to bind cargo to ONE
    Subunit during transport. <site> on all other Subunits remains
    empty.
    """
    assert site in Source.monomer.sites_dict, \
        "Required site %s not present in %s as required"%(site, Source.monomer.name)
    assert site in Dest.monomer.sites_dict, \
        "Required site %s not present in %s as required"%(site, Dest.monomer.name)

    for i in range(min_size, max_size+1):
        # require all ringp subunit <tsite> sites to be empty for Ringp match
        Ringp = ringp_species(Subunit({site: None}), i)

        #r1name = '%s_ringp_%d_transport_%s_cplx' % (Source.monomer.name, i, Subunit.monomer.name)
        r1name = Source.monomer.name + "_" + Subunit.monomer.name + "_".join(["".join(map(str, j)) for j in Subunit.site_conditions.iteritems()]) + "_" + str(i) + "mer_cplx"
        #r2_name = '%s_ringp_%d_transport_%s_dssc' % (Source.monomer.name, i, Subunit.monomer.name)
        r2name = Source.monomer.name + "_" + Subunit.monomer.name + "_".join(["".join(map(str, j)) for j in Subunit.site_conditions.iteritems()]) + "_" + str(i) + "mer_dssc"

        rule_rates = rates[i-min_size]
        CRingp = Ringp.copy()
        tbondnum = i + 1
        CRingp.monomer_patterns[0].site_conditions[site] = tbondnum
        Complex = CRingp % Source({site: tbondnum})
        Rule(r1name, Ringp + Source({site: None}) <> Complex, *rule_rates[0:2])
        Rule(r2name, Complex >> Ringp + Dest({site: None}), rule_rates[2])

def one_step_conv(Sub1, Sub2, Prod, klist, site='bf'):
    """ Convert two Sub species into one Prod species:
    Sub + Sub <> Prod
    """
    kf, kr = klist
    r1_name = 'conv_%s_%s_to_%s'%(Sub1.monomer.name, Sub2.monomer.name, Prod.monomer.name)
    assert site in Sub1.monomer.sites_dict, \
        "Required site %s not present in %s as required"%(site, Sub.monomer.name)
    assert site in Sub2.monomer.sites_dict, \
        "Required site %s not present in %s as required"%(site, Sub.monomer.name)
    # create the sites for the monomers

    Sub1.site_conditions[site] = None
    Sub2.site_conditions[site] = None

    # combine the monomers into a product step rule
    Rule(r1_name, Sub1 + Sub2 <> Prod, kf, kr)
    
def simple_bind(Sub1, Sub2, klist, site='bf'):
    """Automation of the Sub1 + Sub2 <> Sub1:Sub2 one-step complex formation. 
    This function assumes that there is a site named 'bf' which, for simplicity
    need not be passed"""
    
    kf, kr = klist
    
    # FIXME: this will fail if the argument passed is a complex... 
    r1_name = 'cplx_%s_%s' % (Sub1.monomer.name, Sub2.monomer.name)
    
    assert site in Sub1.monomer.sites_dict, \
        "Required site %s not present in %s as required"%(site, Sub1.monomer.name)
    assert site in Sub2.monomer.sites_dict, \
        "Required site %s not present in %s as required"%(site, Sub2.monomer.name)
    
    # create the site conditions for the complex
    s1tmpdict = Sub1.site_conditions.copy()
    s2tmpdict = Sub2.site_conditions.copy()
    
    s1tmpdict[site] = 1
    s2tmpdict[site] = 1

    Sub1Cplx = Sub1.monomer(s1tmpdict)
    Sub2Cplx = Sub2.monomer(s2tmpdict)

    # create the sites for the monomers
    Sub1.site_conditions[site] = None
    Sub2.site_conditions[site] = None
    # now that we have the complex elements formed we can write the first step rule
    Rule(r1_name, Sub1 + Sub2 <> Sub1Cplx % Sub2Cplx, kf, kr)

inhibit = simple_bind #alias for simplebind

#FIXME: pass klist of sorts?
def simple_bind_table(bindtable, parmlist, lmodel, site='bf'):
    """This assumes that the monomers passed are in their desired state without
    the 'bf' site, which will be used for binding.
    bindtable is a list of lists denoting the reactions between two types of reactants
    as follows:

    bindtable[0]: [                  reactypeA0(args), reactypeA1(args)... reactypeAN(args)]
    bindtable[1]: [reactypeB0(args), 'parmfamA0B0',    'parmfamA1B0'...    'parmfamANB0'   ]
    bindtable[2]: [reactypeB1(args), 'parmfamA0B1',    'parmfamA1B1'...    'parmfamANB1'   ]
    
    the variable 'lmodel' is the model passed for local lookup of parameter variables
    """

    # parse the list, extract reactants, products and parameter families
    #first line is one set of reactants
    react0 = bindtable[0]
    react0st = bindtable[1]
    react1 = [row[0] for row in bindtable[2:]]
    react1st = [row[1] for row in bindtable[2:]]

    # Notice this makes intrxns of size/index intrxns[react1][react0]
    intrxns = [row[2:] for row in bindtable[2:]]
    
    # Add the bf sites to the reactant states dict
    # NOTE: this will reset the value if it is already set.
    # Build the prod states dict from react dicts, change bf to 1
    prod0st = []
    prod1st = []
    for d in react0st:
        d[site] = None
        prod0st.append(d.copy())
    for d in react1st:
        d[site] = None
        prod1st.append(d.copy())
    for d in prod0st:
        d[site] = 1
    for d in prod1st:
        d[site] = 1
    
    # loop over interactions
    pc = 0 # parameter counter, cheap way of keeping track of which param set in the list to use
    for i in range(0, len(react1)):
        for j in range(0, len(react0)):
            if intrxns[i][j] is True:
                # get the parameters from the parmlist
                kf = parmlist[pc][0]
                kr = parmlist[pc][1]
                # rule name
                rname = 'cplx_%s_%s' % (react1[i].name, react0[j].name)
                # create the rule
                #print "Generating  %s:%s complex"%(react1[i].name, react0[j].name)
                Rule(rname, react1[i](react1st[i]) + react0[j](react0st[j]) <>
                     react1[i](prod1st[i]) % react0[j](prod0st[j]), 
                     kf, kr)
                pc += 1
    if pc != len(parmlist):
        print "WARNING, unassigned parameters from list", parmlist
        print "Assigned",pc,"parameter pairs from a total of", len(parmlist)
    

#-------------------------------------------------------------------------
# Random little helper funcs that make it easier to interact w the model.
#-------------------------------------------------------------------------

def get_param_num(model, name):
    for i in range(len(model.parameters)):
        if model.parameters[i].name == name:
            print i, model.parameters[i]
            break
    return i

def plotoutput(simout, norm=True):
    """ Assume norm is true for now
    """
    pylab.ion()
    pylab.figure()
    nplots = len(simout.shape[0] -1)
    
    
    for i in range(nplots): #assume simout[0] is time
        pass
        
