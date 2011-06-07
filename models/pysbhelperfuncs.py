from pysb import *
def twostepmod(Enz, Sub, Prod, kf, kr, kc):
    """Automation of the Enz + Sub <> Enz:Sub >> Enz + Prod two-step catalytic reaction.
    This function assumes that there is a site named 'bf' (bind site for fxn)
    which it uses by default. This also assume Enz returns to its original state.
    In an aim for simplicity, site 'bf' need not be passed when calling the function by
    the reactants, but THE FULL STATE OF THE PRODUCT must be specified"""
    
    # FIXME: this will fail if the argument passed is a Complex object. 

    r1_name = 'cplx_%s_%s' % (Sub.monomer.name, Enz.monomer.name)
    r2_name = 'dssc_%s_via_%s' % (Prod.monomer.name, Enz.monomer.name)
    
    assert 'bf' in Enz.monomer.sites_dict, \
        "Required site 'bf' not present in %s as required"%(Enz.monomer.name)
    assert 'bf' in Sub.monomer.sites_dict, \
        "Required site 'bf' not present in %s as required"%(Sub.monomer.name)

    # make the intermediate complex components
    etmpdict = Enz.site_conditions.copy()
    stmpdict = Sub.site_conditions.copy()
    
    etmpdict['bf'] = 1
    stmpdict['bf'] = 1

    EnzCplx = Enz.monomer(etmpdict)
    SubCplx = Sub.monomer(stmpdict)

    # add the 'bf' site to the patterns
    Enz.site_conditions['bf'] = None
    Sub.site_conditions['bf'] = None

    # now that we have the complex elements formed we can write the first step rule
    Rule(r1_name, Enz + Sub <> EnzCplx % SubCplx, kf, kr)
    
    # and finally the rule for the catalytic transformation
    Rule(r2_name, EnzCplx % SubCplx >> Enz + Prod, kc)

def twostepconv(Sub1, Sub2, Prod, kf, kr, kc):
    """Automation of the Sub1 + Sub2 <> Sub1:Sub2 >> Prod two-step reaction (i.e. dimerization).
    This function assumes that there is a site named 'bf' (bind site for fxn)
    which it uses by default. Site 'bf' need not be passed when calling the function."""
    
    # FIXME: this will fail if the argument passed is a Complex

    r1_name = 'cplx_%s_%s' % (Sub2.monomer.name, Sub1.monomer.name)
    r2_name = 'cplx_%s_via_%s__%s' % (Prod.monomer.name, Sub1.monomer.name, Sub2.monomer.name)
    
    assert 'bf' in Sub1.monomer.sites_dict, \
        "Required site 'bf' not present in %s as required"%(Sub1.monomer.name)
    assert 'bf' in Sub2.monomer.sites_dict, \
        "Required site 'bf' not present in %s as required"%(Sub2.monomer.name)

    # make the intermediate complex components
    s1tmpdict = Sub1.site_conditions.copy()
    s2tmpdict = Sub2.site_conditions.copy()
    
    s1tmpdict['bf'] = 1
    s2tmpdict['bf'] = 1

    Sub1Cplx = Sub1.monomer(s1tmpdict)
    Sub2Cplx = Sub2.monomer(s2tmpdict)

    # add the 'bf' site to the patterns
    Sub1.site_conditions['bf'] = None
    Sub2.site_conditions['bf'] = None

    # now that we have the complex elements formed we can write the first step rule
    Rule(r1_name, Sub1 + Sub2 <> Sub1Cplx % Sub2Cplx, kf, kr)
    
    # and finally the rule for the catalytic transformation
    Rule(r2_name, Sub1Cplx % Sub2Cplx >> Prod, kc)

def simpledim(Sub, Prod, kf, kr):
    """ Convert two Sub species into one Prod species:
    Sub + Sub <> Prod
    """
    r1_name = 'dimer_%s_to_%s'%(Sub.monomer.name, Prod.monomer.name)
    assert 'bf' in Sub.monomer.sites_dict, \
        "Required site 'bf' not present in %s as required"%(Sub.monomer.name)

    # create the sites for the monomers
    Sub.site_conditions['bf'] = None

    # combine the monomers into a product step rule
    Rule(r1_name, Sub + Sub <> Prod, kf, kr)

def onestepconv(Sub1, Sub2, Prod, kf, kr):
    """ Convert two Sub species into one Prod species:
    Sub + Sub <> Prod
    """
    r1_name = 'conv_%s_%s_to_%s'%(Sub1.monomer.name, Sub2.monomer.name, Prod.monomer.name)
    assert 'bf' in Sub1.monomer.sites_dict, \
        "Required site 'bf' not present in %s as required"%(Sub.monomer.name)
    assert 'bf' in Sub2.monomer.sites_dict, \
        "Required site 'bf' not present in %s as required"%(Sub.monomer.name)
    # create the sites for the monomers

    Sub1.site_conditions['bf'] = None
    Sub2.site_conditions['bf'] = None

    # combine the monomers into a product step rule
    Rule(r1_name, Sub1 + Sub2 <> Prod, kf, kr)
    
def simplebind(Sub1, Sub2, kf, kr):
    """Automation of the Sub1 + Sub2 <> Sub1:Sub2 one-step complex formation. 
    This function assumes that there is a site named 'bf' which, for simplicity
    need not be passed"""
    
    # FIXME: this will fail if the argument passed is a complex... 
    r1_name = 'cplx_%s_%s' % (Sub1.monomer.name, Sub2.monomer.name)
    
    assert 'bf' in Sub1.monomer.sites_dict, \
        "Required site 'bf' not present in %s as required"%(Sub1.monomer.name)
    assert 'bf' in Sub2.monomer.sites_dict, \
        "Required site 'bf' not present in %s as required"%(Sub2.monomer.name)
    
    # create the site conditions for the complex
    s1tmpdict = Sub1.site_conditions.copy()
    s2tmpdict = Sub2.site_conditions.copy()
    
    s1tmpdict['bf'] = 1
    s2tmpdict['bf'] = 1

    Sub1Cplx = Sub1.monomer(s1tmpdict)
    Sub2Cplx = Sub2.monomer(s2tmpdict)

    # create the sites for the monomers
    Sub1.site_conditions['bf'] = None
    Sub2.site_conditions['bf'] = None
    # now that we have the complex elements formed we can write the first step rule
    Rule(r1_name, Sub1 + Sub2 <> Sub1Cplx % Sub2Cplx, kf, kr)

# def getparm(plist, name):
#     for i in plist:
#         if i.name == name:
#             return i
#         else:
#             continue
#     print "PARAMETER %s NOT FOUND"%(name)
#     LookupError

def sbindtable(bindtable, lmodel):
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
        d['bf'] = None
        prod0st.append(d.copy())
    for d in react1st:
        d['bf'] = None
        prod1st.append(d.copy())
    for d in prod0st:
        d['bf'] = 1
    for d in prod1st:
        d['bf'] = 1
    
    # loop over interactions
    for i in range(0, len(react1)):
        for j in range(0, len(react0)):
            if intrxns[i][j] is True:
                # build the name of the forward/reverse parameters
                sparamf = react1[i].name.lower()+react0[j].name.lower()+'f'
                sparamr = react1[i].name.lower()+react0[j].name.lower()+'r'
                kf = lmodel.parameter(sparamf)
                kr = lmodel.parameter(sparamr)
                # rule name
                rname = 'cplx_%s_%s' % (react1[i].name, react0[j].name)
                # create the rule
                print "Generating  %s:%s complex"%(react1[i].name, react0[j].name)
                Rule(rname, react1[i](react1st[i]) + react0[j](react0st[j]) <>
                     react1[i](prod1st[i]) % react0[j](prod0st[j]), 
                     kf, kr)
    
    
