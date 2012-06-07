import pysb.macros as macros

site_name = 'b'
def declare_monomers():
    Monomer('tBid', ['b'])
    Monomer('Bax', ['b', 'state'], {'state':['C', 'A']})
    Monomer('Bcl2', ['b'])
    Monomer('Bad', ['b'])
    Monomer('Pore')

def bind(a, b, klist):
    return macros.bind(a, site_name, b, site_name, klist)

def bind_table(table):
    return macros.bind(table, row_site, col_site)

def pore_assembly(monomer, size, pore, klist):
    #Rule
    pass
