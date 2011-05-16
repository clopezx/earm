from pysb import *
def catact(Enz, Sub, Prod, kf, kr, kc):
    """Automation of the Enz + Sub <> Enz:Sub >> Enz+Prod two-step catalytic reaction
    this function assumes that there is a site named 'bf' (bind catalysis)
    hich it uses by default. This also assume Enz returns to its original state."""
    
    # FIXME: this will fail if the argument passed is anything more than a MonomerPattern
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

    # now that we have the complex pattern formed we can write the first step rule
    Rule(r1_name, Enz + Sub <> EnzCplx % SubCplx, kf, kr)
    
    # and finally the rule for the catalytic transformation
    Rule(r2_name, EnzCplx % SubCplx >> Enz + Prod, kc)
    

