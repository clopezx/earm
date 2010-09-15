from pysb import *

# FUNCTIONS
# def can_activate(*monomers):
#     for m in monomers:
#         m.sites += ['active']
#         m.sites_dict['active'] = None
#         m.site_states['active'] = ['no', 'yes']
#

def rdata_list(rulename, rule_func, *rows):
    """This allows the processing of data in multiple rows which makes simple
       consecutive reactions easy"""
    counter = 0
    for r in rows:
        rule_func(rulename, *r)
        counter += 1
    #print "# added", counter, "rules"

def rdata_set(rulename, rule_func, firstrow, *restrows):
    """This allows the processing of data rule sets of arbitrary reactions of the form
    dealing with r1, r2 only"""
    counter = 0
    #print '\nfirstrow_len:', len(firstrow), 'restrows_len:', len(restrows)
    for j in range(0,len(restrows)):
        for i in range(0,len(firstrow)):          
            reac1= firstrow[i]
            reac2= restrows[j][0]
            param= restrows[j][i+1]
            #print '##r1:', reac1
            #print '##r2:', reac2
            #print param
            rule_func(rulename, reac1, reac2, param)
            counter += 1
    #print "# added", counter, "rules"

# can_activate(CASP8, BIDC, BIDM, BAXC, BAXM, BAK)
# dimerize(BAXM, k_baxm_dim_f, k_baxm_dim_r)
# dimerize(BAK,  k_bak_dim_f,  k_bak_dim_r)
#
# def activate(monomerI1, monomerI2, monomerO, ratef, rater, ratea):

#def simpleComplex(monomer1, monomer2):
#    """This automates the  A + B <-> A:B rule"""

#def activeComplex(monomer1, monomer2):
#    """This automates the A + B <-> A:B -> [AB]* rule"""
    
#def active_w_enzyme(monomer1, enzyme):
#    """This automates the A + B <-> A:B -> A* + B rule, where B is an enzyme"""
    
Model()

# COMPARTMENTS
# I am placing this here just as a reminder that it needs to be dealt with
#(defcon cell-location (location)
#  (&optional (id := *name*) 
#	     &property
#	     (cytoplasm compartment :#= [[compartment] {.size.value := 1.0}]); 1000.0, um^3 cell volume of 1pL (in um^3)
#	     (dish      compartment :#= [[compartment] {.size.value := 1.0}]); 1.0E6, um^3 1E-9L = dish volume 
#	     (mitochondria compartment :#= [compartment])
#	     (cell-membrane membrane :#= (let ((d .dish)
#					       (c .cytoplasm))
#					   [[membrane] :outer d :inner c {.size.value := 1.0}])); 706 um^2, 7.5um cell radius
#	     (mito-membrane membrane :#= (let ((e .mitochondria)
#					       (c .cytoplasm))
#					   [[membrane] :inner e :outer c {.size.value := 1.0}])); 900 um^2, 300 mito, 3um^2/mito
#	     (inverse-mito-membrane membrane :#= .mito-membrane.inverse))) ; invert inner and outer
#
# ;; DEFINE CELL AND SPECIES IN COMPARTMENTS
# (define cell [cell-location])
# cell.dish.(contains [TRAIL])
# cell.cell-membrane.(contains [DRECEPTOR])
# cell.cytoplasm.(contains [caspase8] [caspase3] [caspase6] [caspase9] [BIDC] [BAXC] [BIM] 
# 			 [BCLXLC]   [MCL1C]    [NOXA]     [BAD]      [MULE] [XIAP] [PARP] 
# 			 [APAF]     [apopsome])
# cell.mito-membrane.(contains [BIDM] [BAXM] [BAK] [BCLXLM] [BCL2] [MCL1M] [BAXMPORE] [BAKPORE])
# cell.mitochondria.(contains [CYTC] [SMAC])
#

# ***CONCENTRATIONS***
# Concentrations in number per cell
# EXAMPLE: Parameter('EGF_tot',  1.2e6)
#
Parameter('TRAIL_0',     3000) # molecules/cell original: 3000
Parameter('DReceptor_0',  200) # molecules/cell
Parameter('Csp8_0',       2E4) # molecules/cell
Parameter('Csp3_0',       1E4) # molecules/cell
Parameter('Csp6_0',       1E4) # molecules/cell
Parameter('Csp9_0',       1E5) # molecules/cell
Parameter('Bid_0',        4E4) # molecules/cell
Parameter('Bax_0',        1E5) # molecules/cell
Parameter('Bim_0',          0) # molecules/cell 5E4
Parameter('Bak_0',        5E4) # molecules/cell
Parameter('BclXl_C_0' , 1.98E4) # molecules/cell
Parameter('BclXl_M_0' ,    200) # molecules/cell
Parameter('Bcl2_0',       2E4) # molecules/cell
Parameter('Mcl1_C_0',   1.98E4) # molecules/cell
Parameter('Mcl1_M_0',      200) # molecules/cell
Parameter('Bad_0',        5E4) # molecules/cell
Parameter('Noxa_0',       5E4) # molecules/cell
Parameter('MULE_0',       1E4) # molecules/cell
Parameter('XIAP_0',       1E5) # molecules/cell
Parameter('Parp_0',       1E6) # molecules/cell
Parameter('Apaf_0',       1E5) # molecules/cell
Parameter('CytC_0',       5E5) # molecules/cell
Parameter('Smac_0',       1E5) # molecules/cell

# ***RATE CONSTANTS***
# Biomolecular rate constants are in #/cell
# Bimolecular rates in (#/cell)^-1s^-1
# Unimolecular rates in s^-1
# EXAMPLE: Parameter('kp1', 1.67e-06)
Parameter('KBIDCCSP8F',	       3.0E-8)   # (#/cell)^-1 sec^-1
Parameter('KBIDCCSP8R',	       1.0E-3)   # sec^-1
Parameter('KBIDCCSP8C',	       1.0   )   # sec^-1
Parameter('KBIDCBAXCF',	       4.5E-8)   # (#/cell)^-1 sec^-1
Parameter('KBIDCBAXCR',        1.0E-3)   # sec^-1
Parameter('KBIDCBAXCC',        1.0   )   # sec^-1
Parameter('KBIDMBAXCF',        2.2E-8)   # (#/cell)^-1 sec^-1
Parameter('KBIDMBAXCR',        1.0E-3)   # sec^-1
Parameter('KBIDMBAXCC',        1.0 )     # sec^-1
Parameter('KBIDCBAKF', 	       2.0E-6 )  # (#/cell)^-1 sec^-1
Parameter('KBIDCBAKR', 	       1.0E-3 )  # sec^-1
Parameter('KBIDCBAKC', 	       1.0 )     # sec^-1
Parameter('KBIDMBAKF', 	       2.0E-6 )  # (#/cell)^-1 sec^-1
Parameter('KBIDMBAKR', 	       1.0E-3 )  # sec^-1
Parameter('KBIDMBAKC', 	       1.0 )     # sec^-1
Parameter('KBIMBAXCF', 	       0.0 )     # (#/cell)^-1 sec^-1 2.4E-7
Parameter('KBIMBAXCR', 	       0.0 )     # sec^-1 	       1.0E-3
Parameter('KBIMBAXCC', 	       0.0 )     # sec^-1 	       1.0
Parameter('KBIMBAKF', 	       0.0 )     # (#/cell)^-1 sec^-1 2.0E-6
Parameter('KBIMBAKR', 	       0.0 )     # sec^-1 	       1.0E-3
Parameter('KBIMBAKC', 	       0.0 )     # sec^-1             1.0
Parameter('KBAXAUTOACTF',      0.0 )     # sec^-1 *****NOT S
Parameter('KBAXAUTOACTR',      0.0 )     # sec^-1 *****NOT S
Parameter('KBAXMDIMF', 	       1.0e-6 )  # (#/cell)^-1 sec^-1
Parameter('KBAXMDIMR', 	       1.0e-3 )  # sec^-1
Parameter('KBAXMTETF', 	       1.0e-6 )  # (#/cell)^-1 sec^-1
Parameter('KBAXMTETR', 	       1.0e-3 )  # sec^-1
Parameter('KBAXMPORE',         1.0 )     # s
Parameter('KBAKAUTOACTF',      0.0 )     # sec^-1 *****NOT S
Parameter('KBAKAUTOACTR',      0.0 )     # sec^-1 *****NOT S
Parameter('KBAKDIMF', 	       2.0e-6 )  # (#/cell)^-1 sec^-1
Parameter('KBAKDIMR', 	       1.0e-3 )  # sec^-1
Parameter('KBAKTETF', 	       2.0e-6 )  # (#/cell)^-1 sec^-1
Parameter('KBAKTETR', 	       1.0e-3 )  # sec^-1
Parameter('KBAKPORE', 	       1.0 )     # sec^-1
Parameter('KBCLXLCBAXCF',      8.3e-10 ) # (#/cell)^-1 sec^-1
Parameter('KBCLXLCBAXCR',      1.0e-3 )  # sec^-1
Parameter('KBCLXLMBAXCF',      8.3e-10 ) # (#/cell)^-1 sec^-1
Parameter('KBCLXLMBAXCR',      1.0e-3 )  # sec^-1
Parameter('KBCLXLCBAXMF',      8.3e-10 ) # (#/cell)^-1 sec^-1
Parameter('KBCLXLCBAXMR',      1.0e-3 )  # sec^-1
Parameter('KBCLXLMBAXMF',      8.3e-10 ) # (#/cell)^-1 sec^-1
Parameter('KBCLXLMBAXMR',      1.0e-3 )  # sec^-1
Parameter('KMCL1CBAXCF',       1.0e-6 )  # ******C
Parameter('KMCL1CBAXCR',       1.0e-3 )  # ******C
Parameter('KMCL1CBAXMF',       1.0e-6 )  # ******C
Parameter('KMCL1CBAXMR',       1.0e-3 )  # ******C
Parameter('KMCL1MBAXCF',       1.0e-6 )  # ******C
Parameter('KMCL1MBAXCR',       1.0e-3 )  # ******C
Parameter('KMCL1MBAXMF',       1.0e-6 )  # ******C
Parameter('KMCL1MBAXMR',       1.0e-3 )  # ******C
Parameter('KBCL2BAXCF',        1.7e-9 )  # (#/cell)^-1 sec^-1
Parameter('KBCL2BAXCR',        1.0e-3 )  # sec^-1
Parameter('KBCL2BAXMF',        1.0E-6 )  # (#/cell)^-1 sec^-1
Parameter('KBCL2BAXMR',        1.0e-3 )  # sec^-1
Parameter('KBCLXLCBAKF',       3.3e-8 )  # (#/cell)^-1 sec^-1
Parameter('KBCLXLCBAKR',       1.0e-3 )  # sec^-1
Parameter('KBCLXLMBAKF',       3.3e-8 )  # (#/cell)^-1 sec^-1
Parameter('KBCLXLMBAKR',       1.0e-3 )  # sec^-1
Parameter('KMCL1CBAKF',        1.7e-6 )  # (#/cell)^-1 sec^-1
Parameter('KMCL1CBAKR',        1.0e-3 )  # sec^-1
Parameter('KMCL1MBAKF',        1.7e-6 )  # (#/cell)^-1 sec^-1
Parameter('KMCL1MBAKR',        1.0e-3 )  # sec^-1
Parameter('KBCL2BAKF', 	       1.7e-6 )  # ******C
Parameter('KBCL2BAKR', 	       1.0e-3 )  # ******C
Parameter('KBCLXLCBIDCF',      3.2e-7 )  # (#/cell)^-1 sec^-1
Parameter('KBCLXLCBIDCR',      5.5e-3 )  # sec^-1
Parameter('KBCLXLMBIDCF',      3.2e-7 )  # (#/cell)^-1 sec^-1
Parameter('KBCLXLMBIDCR',      5.5e-3 )  # sec^-1
Parameter('KBCLXLCBIDMF',      3.2e-7 )  # (#/cell)^-1 sec^-1
Parameter('KBCLXLCBIDMR',      5.5e-3 )  # sec^-1
Parameter('KBCLXLMBIDMF',      3.2e-7 )  # (#/cell)^-1 sec^-1
Parameter('KBCLXLMBIDMR',      5.5e-3 )  # sec^-1
Parameter('KMCL1CBIDCF',       1.7e-6 )  # (#/cell)^-1 sec^-1
Parameter('KMCL1CBIDCR',       1.0e-3 )  # sec^-1
Parameter('KMCL1CBIDMF',       1.7e-6 )  # (#/cell)^-1 sec^-1
Parameter('KMCL1CBIDMR',       1.0e-3 )  # sec^-1
Parameter('KMCL1MBIDCF',       1.7e-6 )  # (#/cell)^-1 sec^-1
Parameter('KMCL1MBIDCR',       1.0e-3 )  # sec^-1
Parameter('KMCL1MBIDMF',       1.7e-6 )  # (#/cell)^-1 sec^-1
Parameter('KMCL1MBIDMR',       1.0e-3 )  # sec^-1
Parameter('KBCL2BIDCF',        1.0e-6 )  # (#/cell)^-1 sec^-1
Parameter('KBCL2BIDCR',        1.0e-3 )  # sec^-1
Parameter('KBCL2BIDMF',        9.1e-7 )  # (#/cell)^-1 sec^-1
Parameter('KBCL2BIDMR',        1.0e-3 )  # sec^-1
Parameter('KBCLXLCBIMF',       5.3e-7 )  # (#/cell)^-1 sec^-1
Parameter('KBCLXLCBIMR',       2.3e-3 )  # sec^-1
Parameter('KBCLXLMBIMF',       5.3e-7 )  # (#/cell)^-1 sec^-1
Parameter('KBCLXLMBIMR',       2.3e-3 )  # sec^-1
Parameter('KBCL2BIMF', 	       5.0e-9 )  # (#/cell)^-1 sec^-1
Parameter('KBCL2BIMR', 	       1.4e-4 )  # sec^-1
Parameter('KMCL1CBIMF',        1.3e-6 )  # (#/cell)^-1 sec^-1
Parameter('KMCL1CBIMR',        2.2e-3 )  # sec^-1
Parameter('KMCL1MBIMF',        1.3e-6 )  # (#/cell)^-1 sec^-1
Parameter('KMCL1MBIMR',        2.2e-3 )  # sec^-1
Parameter('KBCLXLCBADF',       6.8e-7 )  # (#/cell)^-1 sec^-1
Parameter('KBCLXLCBADR',       6.5e-4 )  # sec^-1
Parameter('KBCLXLMBADF',       6.8e-7 )  # (#/cell)^-1 sec^-1
Parameter('KBCLXLMBADR',       6.5e-4 )  # sec^-1
Parameter('KBCL2BADF', 	       1.3e-7 )  # (#/cell)^-1 sec^-1
Parameter('KBCL2BADR', 	       1.0e-3 )  # sec^-1
Parameter('KMCL1CNOXAF',       1.3e-8 )  # (#/cell)^-1 sec^-1
Parameter('KMCL1CNOXAR',       3.7 )     # sec^-1
Parameter('KMCL1CNOXAC',       0.0001 )  # sec^-1 ****NOT S
Parameter('KMCL1MNOXAF',       1.3e-8 )  # (#/cell)^-1 sec^-1
Parameter('KMCL1MNOXAR',       3.7 )     # sec^-1
Parameter('KMCL1MNOXAC',       0.0001 )  # sec^-1 ****NOT S
Parameter('KMCL1CMULEF',       2.0e-6 )  # (#/cell)^-1 sec^-1
Parameter('KMCL1CMULER',       1.0e-3 )  # sec^-1
Parameter('KMCL1CMULEC',       0.0001 )  # sec^-1
Parameter('KMCL1MMULEF',       2.0e-6 )  # (#/cell)^-1 sec^-1
Parameter('KMCL1MMULER',       1.0e-3 )  # sec^-1
Parameter('KMCL1MMULEC',       0.0001 )  # sec^-1
Parameter('KCSP3CSP8F',        8.0e-8 )  # (#/cell)^-1 sec^-1
Parameter('KCSP3CSP8R',        1.0e-3 )  # sec^-1
Parameter('KCSP3CSP8C',        1.0 )     # sec^-1
Parameter('KMCL1CCSP8F',       2.0e-6 )  # (#/cell)^-1 sec^-1
Parameter('KMCL1CCSP8R',       1.0e-3 )  # sec^-1
Parameter('KMCL1CCSP8C',       1.0 )     # sec^-1
Parameter('KMCL1MCSP8F',       2.0e-6 )  # (#/cell)^-1 sec^-1
Parameter('KMCL1MCSP8R',       1.0e-3 )  # sec^-1
Parameter('KMCL1MCSP8C',       1.0 )     # sec^-1
Parameter('KMCL1CCSP3F',       2.0e-6 )  # (#/cell)^-1 sec^-1
Parameter('KMCL1CCSP3R',       1.0e-3 )  # sec^-1
Parameter('KMCL1CCSP3C',       1.0 )     # sec^-1
Parameter('KMCL1MCSP3F',       2.0e-6 )  # (#/cell)^-1 sec^-1
Parameter('KMCL1MCSP3R',       1.0e-3 )  # sec^-1
Parameter('KMCL1MCSP3C',       1.0 )     # sec^-1
Parameter('KCYTCBAXPF',        2.0e-6 )  # (#/cell)^-1 sec^-1
Parameter('KCYTCBAXPR',        1.0e-3 )  # sec^-1
Parameter('KCYTCBAXPC',       10.0 )     # sec^-1
Parameter('KSMACBAXPF',        2.0e-6 )  # (#/cell)^-1 sec^-1
Parameter('KSMACBAXPR',        1.0e-3 )  # sec^-1
Parameter('KSMACBAXPC',       10.0 )     # sec^-1
Parameter('KCYTCBAKPF',        2.0e-6 )  # (#/cell)^-1 sec^-1
Parameter('KCYTCBAKPR',        1.0e-3 )  # sec^-1
Parameter('KCYTCBAKPC',       10.0 )     # sec^-1
Parameter('KSMACBAKPF',        2.0e-6 )  # (#/cell)^-1 sec^-1
Parameter('KSMACBAKPR',        1.0e-3 )  # sec^-1
Parameter('KSMACBAKPC',       10.0 )     # sec^-1
Parameter('KBIDCBIDMF',        1.0e2 )   # sec^-1
Parameter('KBIDCBIDMR',        0.1 )     # sec^-1
Parameter('KBAXCBAXMF',        1.0e-2 )  # sec^-1
Parameter('KBAXCBAXMR',        1.0e-2 )  # sec^-1
Parameter('KBCLXLCBCLXLMF',    1.0e-2 )  # sec^-1
Parameter('KBCLXLCBCLXLMR',    1.0e-3 )  # sec^-1
Parameter('KMCL1CMCL1MF',      1.0e-2 )  # sec^-1
Parameter('KMCL1CMCL1MR',      1.0e-3 )  # sec^-1
Parameter('KCSP8ICSP8A',       7.7E-3 )  # sec^-1 *****MADE UP assuming exponential growth -- E
Parameter('KMCL1CAMCL1CD',  1664.0 )     # sec^-1 *****MADE UP assuming exponential decay with half life of 40 m
Parameter('KMCL1MAMCL1MD',  1664.0 )     # sec^-1
Parameter('KTRAILRECF',        4.0e-7 )  # (#/cell)^-1 sec^-1
Parameter('KTRAILRECR',        1.0e-3 )  # sec^-1
Parameter('KTRAILRECC',        1.0e-5 )  # sec^-1
Parameter('KDRECCSP8F',        6.0E-8 )  # (#/cell)^-1 sec^-1
Parameter('KDRECCSP8R',        1.0e-3 )  # sec^-1
Parameter('KDRECCSP8C',        1.0 )     # sec^-1
Parameter('KCSP3CSP6F',        1.0e-6 )  # (#/cell)^-1 sec^-1
Parameter('KCSP3CSP6R',        1.0e-3 )  # sec^-1
Parameter('KCSP3CSP6C',        1.0 )     # sec^-1
Parameter('KCSP3XIAPF',        5.0e-6 )  # (#/cell)^-1 sec^-1
Parameter('KCSP3XIAPR',        1.0e-3 )  # sec^-1
Parameter('KCSP3XIAPC',        0.1 )     # sec^-1
Parameter('KCSP3PARPF',        2.0E-7 )  # (#/cell)^-1 sec^-1
Parameter('KCSP3PARPR',        1.0e-2 )  # sec^-1
Parameter('KCSP3PARPC',        1.0 )     # sec^-1
Parameter('KCSP6CSP8F',        3.0e-8 )  # (#/cell)^-1 sec^-1
Parameter('KCSP6CSP8R',        1.0e-3 )  # sec^-1
Parameter('KCSP6CSP8C',        1.0 )     # sec^-1
Parameter('KCYTCAPAFF',        5.0e-7 )  # sec^-1
Parameter('KCYTCAPAFR',        1.0e-3 )  # sec^-1
Parameter('KCYTCAPAFC',        1.0 )     # sec^-1
Parameter('KAPAFCSP9F',        5.0e-8 )  # sec^-1
Parameter('KAPAFCSP9R',        1.0e-3 )  # sec^-1
Parameter('KAPOPCSP3F',        1.0e-6 )  # sec^-1
Parameter('KAPOPCSP3R',        1.0e-3 )  # sec^-1
Parameter('KAPOPCSP3C',        1.0 )     # sec^-1
Parameter('KAPOPXIAPF',        1.0e-4 )  # sec^-1
Parameter('KAPOPXIAPR',        1.0e-3 )  # sec^-1
Parameter('KSMACXIAPF',        7.0e-6 )  # sec^-1
Parameter('KSMACXIAPR',        1.0e-3 )  # sec^-1

# ***MONOMERS***
# example:
# EGFR(l, r, Y1068~U~P, Y1148~U~P)
# Monomer('EGFR',
#        ['l','r','Y1068','Y1148'],
#        { 'Y1068': ['U','P'],
#          'Y1148': ['U','P'] }
#        )
#
#
Monomer('TRAIL',     ['b']) # trail ligand
Monomer('DReceptor', ['bt', 'b', 'state'], {'state': ['I', 'A']}) # death receptor, compartment: membrane
Monomer('Csp8',      ['b', 'state'],  {'state': ['I', 'A']}) # "caspase 8 with on/off site"
Monomer('Bid',       ['b', 'state', 'loc'], {'state':['I','A'], 'loc':['C','M']}) # "BID, inactive = untruncated"
Monomer('Bax',       ['b', 't1', 't2', 'state', 'loc'], {'state': ['I','A'], 'loc':['C','M']}) # cytoplasmic bax
Monomer('Bak',       ['b', 't1', 't2', 'state'], {'state': ['I','A']}) #compartment: membrane
Monomer('BaxPore',   ['b']) # compartment: membrane
Monomer('BakPore',   ['b']) # compartment: membrane
Monomer('Mcl1',      ['b', 'loc', 'deg'], {'loc':['C','M'], 'deg':['U','D']}) # can degrade, undeg by default
Monomer('BclXl',     ['b', 'loc'], {'loc':['C','M']})
Monomer('Bim',       ['b'])
Monomer('Bcl2',      ['b', 'loc'], {'loc':'M'}) # compartment: membrane ****
Monomer('Csp3',      ['b', 'state'], {'state':['I','A','U']}) #state Inactive, Active, Ubiquinated
Monomer('Bad',        'b')
Monomer('NOXA',       'b')
Monomer('MULE',      ['b', 'state', 'usite'],
                     {'state':['I', 'A'], # default inactive
                      'usite': ['N', 'U']} # non-ubiquitinated/ubiquinated, default non-ubiquinated
        )
Monomer('CASPASE6',  ['CCSP3', 'CCSP8', 'state'], {'state': ['I', 'A']})
Monomer('CYTC',      ['CPORE', 'CAPAF'])
Monomer('SMAC',      ['CPORE', 'CXIAP'])
Monomer('XIAP',      ['CCSP3', 'CAPOP', 'CSMAC'])
Monomer('PARP',      ['CCSP3', 'state'], {'state': ['I', 'A'] })
Monomer('APAF',      ['CCYTC', 'CCSP9', 'state'], {'state': ['I', 'A'] })
Monomer('CASPASE9',   'CAPAF')
Monomer('APOPSOME',  ['CCSP3', 'CXIAP'])


# ***RULEZ***
#
# Example:
# Rule('grb2_bind_sos',
#      Grb2(SH2=None, SH3=None) + Sos(PR=None) <>
#      Grb2(SH2=None, SH3=1)    % Sos(PR=1),
#      kp5, km5)
#


# TRAIL binding to Receptor
Rule('TRAIL_bind_DRECEPTOR',
     TRAIL(b=None) + DReceptor(bt=None, b=None, state='I') <>
     TRAIL(b=1)    % DReceptor(bt=1,    b=None, state='I'),
     KTRAILRECF, KTRAILRECR)
Rule('TRAILDReceptor_act',
     TRAIL(b=1)    % DReceptor(bt=1, b=None, state='I') >>
     TRAIL(b=1)    % DReceptor(bt=1, b=None, state='A'),  # NOTICE TRAIL STAYS BOUND AT bt!!!
     KTRAILRECC)

# ---->***FLIP INHIBITION NEEDS TO BE ADDED***

# CSP 8 activation
Rule('Csp8_bind_DReceptor',
     DReceptor(b=None, state='A') + Csp8(b=None, state='I') <> # Ignoring Trail bound, implies it is there for Drec to be active
     DReceptor(b=1, state='A')    % Csp8(b=1, state='I'),
     KDRECCSP8F, KDRECCSP8R)
Rule('Csp8_activation',
     DReceptor(b=1, state='A')    % Csp8(b=1, state='I') >>
     DReceptor(b=None, state='A') + Csp8(b=None, state='A'), # cytosolic C8 is released from the 'brane
     KDRECCSP8C)
     
rdata_set("C8_block",
          lambda rname, r1, r2, param:
              (Rule(rname + '_%s_%s_bind' % (r1.name, r2[0].name),
                    r1(b=None, state='A') + r2[0](r2[1], b=None, state='I') <>
                    r1(b=1,    state='A') % r2[0](r2[1], b=1,    state='I'),
                    param[0], param[1]),
               Rule(rname + '_%s_%s_act' % (r1.name, r2[0].name),
                    r1(b=1, state='A') % r2[0](r2[1], b=1, state='I') >>
                    r1(b=None, state='A') + r2[0](r2[1], b=None, state='A'),
                    param[2])),
          [                     Csp8                            ], #r1
          [(Bid, {'loc':'C'}), [KBIDCCSP8F, KBIDCCSP8R, KBIDCCSP8C]],
          [(Csp3, None),       [KCSP3CSP8F, KCSP3CSP8R, KCSP3CSP8C]]
          ) #r2


# BH3 activators binding bax and bak FIXME: CHANGE TO REAL COMPARTMENTS
rdata_set("BH3_activators",
          lambda rname, r1, r2, param:
              (Rule(rname + '_%s_%s%s' % (r1[0].name, r2[0].name, r2[1].get('loc','')),
                    r1[0](r1[1], t1=None, t2=None, state='I', b=None)  + r2[0](r2[1], b=None)   <>
                    r1[0](r1[1], t1=None, t2=None, state='I', b=1   )  % r2[0](r2[1], b=1   ), 
                    param[0], param[1]),
               Rule(rname + '_%s_%s%s_act' % (r1[0].name, r2[0].name, r2[1].get('loc','')),
                    r1[0](r1[1], t1=None, t2=None, state='I', b=1   )  % r2[0](r2[1], b=1   )   >>
                    r1[0](r1[1], t1=None, t2=None, state='A', b=None)  + r2[0](r2[1], b=None),
                    param[2])
               ),
          [                                 (Bax, {'loc':'C'}),                   (Bak, None)                     ], #r1
          [(Bid, {'state':'A', 'loc':'C'}), [KBIDCBAXCF, KBIDCBAXCR, KBIDCBAXCC], [KBIDCBAKF, KBIDCBAKR, KBIDCBAKC]],
          [(Bid, {'state':'A', 'loc':'M'}), [KBIDMBAXCF, KBIDMBAXCR, KBIDMBAXCC], [KBIDMBAKF, KBIDMBAKR, KBIDMBAKC]], 
          [(Bim, {}),                       [KBIMBAXCF,  KBIMBAXCR,  KBIMBAXCC ], [KBIMBAKF,  KBIMBAKR,  KBIMBAKC ]],
          #r2
          )



# Translocations FIXME: CHANGE TO REAL COMPARTMENTS
rdata_list("translocate",
           lambda rname, r1, r2, data:
               Rule(rname + '_%s%s_%s%s' % (r1[0].name, r1[1].get('loc',''), r2[0].name, r2[1].get('loc','')),
                    r1[0](r1[1], b=None) >> r2[0](r2[1], b=None),
                    data[0], data[1]),
           [(Bid,    {'state':'A', 'loc':'C'}), (Bid,   {'state':'A', 'loc':'M'}), [KBIDCBIDMF,     KBIDCBIDMR]],
           [(Bax,    {'state':'A', 'loc':'C'}), (Bax,   {'state':'A', 'loc':'M'}), [KBAXCBAXMF,     KBAXCBAXMR]],
           [(BclXl,  {'loc':'C'}),              (BclXl, {'loc':'M'}), [KBCLXLCBCLXLMF, KBCLXLCBCLXLMR]],
           [(Mcl1,   {'loc':'C'}),              (Mcl1,  {'loc':'M'}), [KMCL1CMCL1MF,   KMCL1CMCL1MR]]
           ) 

#
# AUTOACTIVATION SKIPPED
#
# ;; BAK/BAX "spontaneously" activates itself (australian model)
# (with-substitution-table
#  (($R1   $Kf            $Kr                         )
#   (BAXC  KBAXAUTOACTF   KBAXAUTOACTR                )
#   (BAK   KBAKAUTOACTF   KBAKAUTOACTR                ))
#   [{[$R1 state.inactive __] <<->> [$R1 state.active __]}
#  (.set-rate-function 'mass-action :fwd  $Kf :rev $Kr)])


# OLIGOMERIZATION NOTE:
# Change this format to a single BAX_DIM and BAX_TETRA monomer with multiple
# sites to bind BCL2 later on. For now leave like this to have a 1:1 correspondence
# with the original little-b code
#
# BAK/BAX Oligomerization (dimer formation) FIXME: THIS TAKES PLACE IN THE RIGHT COMPARTMENT
rdata_list("dimer_formation",
           lambda rname, r1, data:
               Rule(rname + '_' + r1[0].name,
                    r1[0](r1[1], t1=None, t2=None, b=None) + r1[0](r1[1], t1=None, t2=None, b=None) <> 
                    r1[0](r1[1], t1=1,    t2=None, b=None) % r1[0](r1[1], t1=None, t2=1,    b=None),
                    data[0], data[1]),
           [(Bax, {'state':'A', 'loc':'M'}), [KBAXMDIMF,  KBAXMDIMR]],
           [(Bak, {'state':'A'}),            [KBAKDIMF,   KBAKDIMR]]
           )

# BAX/BAXM Oligomerization (tetramer formation)
rdata_list("tetramer_formation",
           lambda rname, r1, data:
               Rule(rname + '_' + r1[0].name,
                    r1[0](r1[1], t1=1,    t2=None, b=None)    % r1[0](r1[1], t1=None, t2=1,    b=None) +
                    r1[0](r1[1], t1=2,    t2=None, b=None)    % r1[0](r1[1], t1=None, t2=2,    b=None) <>
                    r1[0](r1[1], t1=1,    t2=3,    b=None)    % r1[0](r1[1], t1=4,    t2=1,    b=None) %
                    r1[0](r1[1], t1=2,    t2=4,    b=None)    % r1[0](r1[1], t1=3,    t2=2,    b=None),
                    data[0], data[1]),
           [(Bax, {'state':'A', 'loc':'M'}), [KBAXMTETF,  KBAXMTETR]],
           [(Bak, {'state':'A'}),        [KBAXMTETF,  KBAXMTETR]]
           )

# BAK/BAXM tetramer becomes a "pore"
rdata_list("tetramer_to_pore",
           lambda rname, r1, r2, data:
               Rule(rname + '_' + r1[0].name,
                    r1[0](r1[1], t1=1,    t2=3,    b=None)    % r1[0](r1[1], t1=4,    t2=1,    b=None) %
                    r1[0](r1[1], t1=2,    t2=4,    b=None)    % r1[0](r1[1], t1=3,    t2=2,    b=None) >>
                    r2(b=None),
                    data[0], data[1]
                    ),
           [(Bax, {'state':'A', 'loc':'M'}), BaxPore,  [KBAXMPORE,   KBAXMPORE]],
           [(Bak, {'state':'A'}),            BakPore,  [KBAKPORE,    KBAKPORE]]
           )

### ANTIAPOPTOTIC inhibiting BAX, BAK FIXME: THIS GOES INTO THE PROPER COMPARTMENTS
rdata_set("anti_apop_bind",
          lambda rname, r1, r2, data:
              Rule(rname + '_%s%s_%s%s' % (r1[0].name, r1[1].get('loc',''), r2[0].name, r2[1].get('loc','')),
                   r1[0](r1[1], b=None) + r2[0](r2[1], b=None) <>
                   r1[0](r1[1], b=1)    % r2[0](r2[1], b=1), 
                   data[0], data[1]
                   ),
          [                        (Bax, {'state':'A', 'loc':'C'}), (Bax, {'state':'A', 'loc':'M'}), (Bak, {'state':'A'})], #r1
          [(BclXl, {'loc':'C'}),   [KBCLXLCBAXCF, KBCLXLCBAXCR], 
                                   [KBCLXLCBAXMF, KBCLXLCBAXMR],
                                   [KBCLXLCBAKF,  KBCLXLCBAKR]],
          [(BclXl, {'loc':'M'}),   [KBCLXLMBAXCF, KBCLXLMBAXCR],      
                                   [KBCLXLMBAXMF, KBCLXLMBAXMF],      
                                   [KBCLXLMBAKF,  KBCLXLMBAKR]],
          [(Bcl2,  {}),            [KBCL2BAXCF,   KBCL2BAXCR],      
                                   [KBCL2BAXMF,   KBCL2BAXMF],      
                                   [KBCL2BAKF,    KBCL2BAKR]],
          [(Mcl1, {'loc':'C'}),    [KMCL1CBAXCF,  KMCL1CBAXCR],      
                                   [KMCL1CBAXMF,  KMCL1CBAXMF],      
                                   [KMCL1CBAKF,   KMCL1CBAKR]],
          [(Mcl1, {'loc':'M'}),    [KMCL1MBAXCF,  KMCL1MBAXCR],      
                                   [KMCL1MBAXMF,  KMCL1MBAXMF],      
                                   [KMCL1MBAKF,   KMCL1MBAKR]]
          )#r2

### ANTIAPOPTOTIC binding BH3s
def helper_anti_apop_bh3s(rname, r1, r2, data):
    if data:
        if isinstance(r1, Monomer):
            r1 = (r1, {})
        print "r1=", r1
        print "r2=", r2
        print "--------------------------------------------"
        Rule(rname + '%s%s_%s%s' % (r1[0].name, r1[1].get('loc',''), r2[0].name, r2[1].get('loc','')),
             r1[0](r1[1], b=None) + r2[0](r2[1], b=None) <>
             r1[0](r1[1], b=1)    % r2[0](r2[1], b=1), 
             data[0], data[1])

rdata_set("anti_apop_bh3s",
          helper_anti_apop_bh3s,
          [                        (Bid, {'state':'A', 'loc':'C'}), (Bid, {'state':'A', 'loc':'M'}), Bim, Bad, NOXA, MULE], #r1
          [(BclXl, {'loc':'C'}),  [KBCLXLCBIDCF, KBCLXLCBIDCR],  # with BIDC
                                  [KBCLXLCBIDMF, KBCLXLCBIDMR],  # with BIDM
                                  [KBCLXLCBIMF,  KBCLXLCBIMR],   # with BIM
                                  [KBCLXLCBADF,  KBCLXLCBADR],   # with BAD
                                  None,                            # with NOXA
                                  None],                           # with MULE
          [(BclXl, {'loc':'M'}),  [KBCLXLMBIDCF, KBCLXLMBIDCR],  # with BIDC
                                  [KBCLXLMBIDMF, KBCLXLMBIDMR],  # with BIDM
                                  [KBCLXLMBIMF,  KBCLXLMBIMR],   # with BIM
                                  [KBCLXLMBADF,  KBCLXLMBADR],   # with BAD
                                  None,                            # with NOXA
                                  None],                           # with MULE
          [(Bcl2, {}),            [KBCL2BIDCF,   KBCL2BIDCR],    # with BIDC
                                  [KBCL2BIDMF,   KBCL2BIDMR],    # with BIDM
                                  [KBCL2BIMF,    KBCL2BIMR],     # with BIM
                                  [KBCL2BADF,    KBCL2BADR],     # with BAD
                                  None,                            # with NOXA
                                  None],                           # with MULE
          [(Mcl1, {'loc':'C'}),   [KMCL1CBIDCF,  KMCL1CBIDCR],   # with BIDC
                                  [KMCL1CBIDMF,  KMCL1CBIDMR],   # with BIDM
                                  [KMCL1CBIMF,   KMCL1CBIMR],    # with BIM
                                  None,                            # with BAD
                                  [KMCL1CNOXAF,  KMCL1CNOXAR],   # with NOXA
                                  [KMCL1CMULEF,  KMCL1CMULER]],  # with MULE
          [(Mcl1, {'loc':'M'}),   [KMCL1MBIDCF,  KMCL1MBIDCR],   # with BIDC
                                  [KMCL1MBIDMF,  KMCL1MBIDMR],   # with BIDM
                                  [KMCL1MBIMF,   KMCL1MBIMR],    # with BIM
                                  None,                            # with BAD
                                  [KMCL1MNOXAF,  KMCL1MNOXAR],   # with NOXA
                                  [KMCL1MMULEF,  KMCL1MMULER]]   # with MULE
          ) #r2

### MCL1 degradation by NOXA/MULE
rdata_set("anti_apop_bind",
          lambda rname, r1, r2, data:
              Rule(rname + '_%s%s_%s%s' % (r1[0].name, r1[1].get('loc',''), r2[0].name, r2[1].get('loc','')),
                   r1[0](r1[1], b=None) + r2[0](r2[1], b=None) <>
                   r1[0](r1[1], b=1)    % r2[0](r2[1], b=1), 
                   data[0], data[1]
                   ),
          [                        (NOXA, {}), (MULE, {})], #r1
          [(Mcl1, {'loc':'C'}),   [], 
                                   []]
          [(Mcl1, {'loc':'C'}),   [], 
                                   []]
          ) #r2
          


### CSP8 activates BidC and CSP3

### CSP3 activates CSP6

### CSP3 gets tagged for ubiquination by XIAP

### CSP3 cleaves PARP and inactivates it

### CSP6 activates CSP8

### CSP8/CSP3 induce MCL1 degradation

### CYTC and SMAC released via BAX/BAK pores

### CYTC complexes with APAF

### APAF and CSP9 form APOP

### APOP binds to XIAP

### SMAC binds to XIAP

#-------------------------------------------------------------------------------
#OBSERVABLES
Observe('tBid', Bid(b=None, state='A', loc='C'))

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
Initial(BclXl(b=None, loc='C'), BclXl_C_0)
Initial(BclXl(b=None, loc='M'), BclXl_M_0)
Initial(Mcl1(b=None, loc='C', deg='U'), Mcl1_C_0)
Initial(Mcl1(b=None, loc='M', deg='U'), Mcl1_M_0)

###
if __name__ == '__main__':
    from pysb.generator.bng import BngGenerator
    gen = BngGenerator(model)
    print gen.get_content()
    print ""
    print "begin actions"
    print "  generate_network({overwrite=>1});"
    print "  simulate_ode({t_end=>21600,n_steps=>360});" # 6 hours, 1-minute steps
    print "end actions"
