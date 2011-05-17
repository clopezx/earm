from pysb import *
def twostepact(Enz, Sub, Prod, kf, kr, kc):
    """Automation of the Enz + Sub <> Enz:Sub >> Enz+Prod two-step catalytic reaction
    this function assumes that there is a site named 'bf' (bind site for fxn)
    which it uses by default. This also assume Enz returns to its original state.
    In an aim for simplicity, site 'bf' need not be passed when calling the function."""
    
    # FIXME: this will fail if the argument passed is a Complex

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

    
