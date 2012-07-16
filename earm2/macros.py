"""This function declares a number of macros that are used by many of the
EARM 2 models. They can be divided into the following four categories:

Macros that are specific to the models in EARM 2, but are used by all of
them. The only macro of this type is

- :py:func:`all_observables`

Aliases to generalized macros in pysb.macros that provide default values
for site names or other arguments. Macros of this type include:

- :py:func:`catalyze`
- :py:func:`bind`
- :py:func:`bind_table`
- :py:func:`assemble_pore_sequential`
- :py:func:`pore_transport`

Macros for mechanisms that appear within the models previously published
by the research group of Pingping Shen (or the model from Howells et al.
(2011), which is derived from one of Shen's models):

- :py:func:`assemble_pore_spontaneous`
- :py:func:`displace`
- :py:func:`displace_reversibly`

Macros for mechanisms that appear within the models described in our
group's earlier work, specifically the models described in Albeck
et al. (2008) PLoS Biology:

- :py:func:`catalyze_convert`
- :py:func:`one_step_conv`
- :py:func:`pore_bind`

----------
"""

import pysb.macros as macros
from pysb import *
from pysb import MonomerPattern, ComplexPattern, ComponentSet
from pysb.util import alias_model_components
import functools

# The default site name to be used for binding reactions
site_name = 'bf'

## Observables declarations ===================
def shared_observables():
    """Declare observables commonly used for the TRAIL pathway.

    Declares truncated (and mitochondrial) Bid, cytosolic (i.e., released)
    Smac, and cleaved PARP.
    """

    alias_model_components()
    # Observables
    # ===========
    # Fig 4B from Albeck observes these, normalizes and inverts them
    # Observe('Bid',   Bid(bf=None, state='U'))
    # Observe('PARP',  PARP(bf=None, state='U'))
    # Observe('Smac',  Smac(bf=None, state='mito'))
    Observable('tBid_',  Bid(state='T'))
    Observable('aSmac_', Smac(state='A'))
    Observable('cPARP_', PARP(state='C'))

## Aliases to pysb.macros =====================
def catalyze(enz, sub, product, klist):
    """Alias for pysb.macros.catalyze with default binding site 'bf'.
    """

    return macros.catalyze(enz, site_name, sub, site_name, product, klist)

def bind(a, b, klist):
    """Alias for pysb.macros.bind with default binding site 'bf'.
    """

    return macros.bind(a, site_name, b, site_name, klist)

def bind_table(table):
    """Alias for pysb.macros.bind_table with default binding sites 'bf'.
    """

    return macros.bind_table(table, site_name, site_name)

def assemble_pore_sequential(subunit, size, klist):
    """Alias for pysb.macros.assemble_pore_sequential with default sites.

    Uses 's1' and 's2' as the sites for subunit-subunit binding in the pore.
    """

    return macros.assemble_pore_sequential(subunit, 's1', 's2', size, klist)

def pore_transport(subunit, csource, cdest, ktable):
    """Alias for pysb.macros.pore_transport with default arguments.

    - Uses 's1' and 's2' as the sites for subunit-subunit binding in the pore
    - Uses 'bf' for the binding site on the pore subunits used to bind cargo
    - Uses a default value of 4 for the size of transport-competent pores
    - Uses 'bf' for the binding site on the cargo used to bind the pore
    """

    return macros.pore_transport(subunit, 's1', 's2', 'bf', 4, 4,
                                csource, 'bf', cdest, ktable)

## Macros for the Shen models
def assemble_pore_spontaneous(subunit, klist):
    """Generate the order-4 assembly reaction 4*Subunit <> Pore.
    """

    def pore_rule_name(rule_expression):
        react_p = rule_expression.reactant_pattern
        mp = react_p.complex_patterns[0].monomer_patterns[0]
        subunit_name = macros._monomer_pattern_label(mp)
        pore_name = mp.monomer.name
        return '%s_to_%s%d' % (subunit_name, mp.monomer.name, 4)

    free_subunit = subunit(s1=None, s2=None)
    macros._macro_rule('spontaneous_pore',
        free_subunit + free_subunit + free_subunit + free_subunit <>
        subunit(s1=1, s2=4) % subunit(s1=2, s2=1) % \
        subunit(s1=3, s2=2) % subunit(s1=4, s2=3),
        klist, ['kf', 'kr'], name_func=pore_rule_name)

def displace(lig1, lig2, target, k):
    """Generate unidirectional displacement reaction L1 + L2:T >> L1:T + L2.

    The signature can be remembered with the following formula:
    "lig1 displaces lig2 from target."
    """

    return macros._macro_rule('displace',
         lig1({site_name:None}) + lig2({site_name:1}) % target({site_name:1}) >>
         lig1({site_name:1}) % target({site_name:1}) + lig2({site_name:None}),
         [k], ['k'])

def displace_reversibly(lig1, lig2, target, klist):
    """Generate reversible displacement reaction L1 + L2:T <> L1:T + L2.

    The signature can be remembered with the following formula:
    "lig1 displaces lig2 from target." The first rate given in
    in klist specifies the forward rate of this reaction; the second
    specifies the reverse rate.
    """

    return macros._macro_rule('displace',
         lig1({site_name:None}) + lig2({site_name:1}) % target({site_name:1}) <>
         lig1({site_name:1}) % target({site_name:1}) + lig2({site_name:None}),
         klist, ['kf', 'kr'])

## Macros for the Albeck models
def catalyze_convert(sub1, sub2, product, klist, site=site_name):
    """Automation of the Sub1 + Sub2 <> Sub1:Sub2 >> Prod two-step reaction.

    Because product is created by the function, it must be fully specified.
    """

    macros._verify_sites(sub1, site)
    macros._verify_sites(sub2, site)

    components = macros._macro_rule('bind',
                             sub1({site: None}) + sub2({site: None}) <>
                             sub1({site: 1}) % sub2({site: 1}),
                             klist[0:2], ['kf', 'kr'])
    components |= macros._macro_rule('convert',
                              sub1({site: 1}) % sub2({site: 1}) >> product,
                              [klist[2]], ['kc'])
    return components

def one_step_conv(sub1, sub2, product, klist, site=site_name):
    """ Bind sub1 and sub2 to form one product: sub1 + sub2 <> product.
    """
    kf, kr = klist

    macros._verify_sites(sub1, site)
    macros._verify_sites(sub2, site)

    return macros._macro_rule('convert',
                       sub1({site: None}) + sub2({site: None}) <> product,
                       klist, ['kf', 'kr'])


def pore_bind(subunit, sp_site1, sp_site2, sc_site, size, cargo, c_site,
              klist):
    """Generate rules to bind a monomer to a circular homomeric pore.

    The pore structure is defined by the `pore_species` macro -- `subunit`
    monomers bind to each other from `sp_site1` to `sp_site2` to form a closed
    ring. The binding reaction takes the form pore + cargo <> pore:cargo.

    Parameters
    ----------
    subunit : Monomer or MonomerPattern
        Subunit of which the pore is composed.
    sp_site1, sp_site2 : string
        Names of the sites where one copy of `subunit` binds to the next.
    sc_site : string
        Name of the site on `subunit` where it binds to the cargo `cargo`.
    size : integer
        Number of subunits in the pore at which binding will occur.
    cargo : Monomer or MonomerPattern
        Cargo that binds to the pore complex.
    c_site : string
        Name of the site on `cargo` where it binds to `subunit`.
    klist : list of Parameters or numbers
        List containing forward and reverse rate constants for the binding
        reaction (in that order). Rate constants should either be both Parameter
        objects or both numbers. If Parameters are passed, they will be used
        directly in the generated Rules. If numbers are passed, Parameters
        will be created with automatically generated names based on <TODO>
        and these parameters will be included at the end of the returned
        component list.
    """

    macros._verify_sites(subunit, sc_site)
    macros._verify_sites(cargo, c_site)

    def pore_bind_rule_name(rule_expression, size):
        # Get ReactionPatterns
        react_p = rule_expression.reactant_pattern
        prod_p = rule_expression.product_pattern
        # Build the label components
        # Pore is always first complex of LHS due to how we build the rules
        subunit = react_p.complex_patterns[0].monomer_patterns[0].monomer
        if len(react_p.complex_patterns) == 2:
            # This is the complexation reaction
            cargo = react_p.complex_patterns[1].monomer_patterns[0]
        else:
            # This is the dissociation reaction
            cargo = prod_p.complex_patterns[1].monomer_patterns[0]
        return '%s_%d_%s' % (subunit.name, size,
                             macros._monomer_pattern_label(cargo))

    components = ComponentSet()
    # Set up some aliases that are invariant with pore size
    subunit_free = subunit({sc_site: None})
    cargo_free = cargo({c_site: None})

    #for size, klist in zip(range(min_size, max_size + 1), ktable):

    # More aliases which do depend on pore size
    pore_free = macros.pore_species(subunit_free, sp_site1, sp_site2, size)

    # This one is a bit tricky. The pore:cargo complex must only introduce
    # one additional bond even though there are multiple subunits in the
    # pore. We create partial patterns for bound pore and cargo, using a
    # bond number that is high enough not to conflict with the bonds within
    # the pore ring itself.
    # Start by copying pore_free, which has all cargo binding sites empty
    pore_bound = pore_free.copy()
    # Get the next bond number not yet used in the pore structure itself
    cargo_bond_num = size + 1
    # Assign that bond to the first subunit in the pore
    pore_bound.monomer_patterns[0].site_conditions[sc_site] = cargo_bond_num
    # Create a cargo source pattern with that same bond
    cargo_bound = cargo({c_site: cargo_bond_num})
    # Finally we can define the complex trivially; the bond numbers are
    # already present in the patterns
    pc_complex = pore_bound % cargo_bound

    # Create the rules
    name_func = functools.partial(pore_bind_rule_name, size=size)
    components |= macros._macro_rule('pore_bind',
                              pore_free + cargo_free <> pc_complex,
                              klist[0:2], ['kf', 'kr'],
                              name_func=name_func)

    return components
