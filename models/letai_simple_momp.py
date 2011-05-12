from pysb import *

# Show how to incorporate complex rules easily
# Based on Letai data

Model()

# Species monomer sections
Monomer('tBid', ['b', 'state'], {'state':['A', 'I']})
Monomer('Bax', ['b', 'm', 'state'], {'state':['A', 'I']})
Monomer('CytC', ['b', 'state'], {'state':['mito', 'cyto']})
Monomer('BaxPore')

# Parameter section
Parameter('[NAME]', [VALUE])

# Bax binds to active tBid, Bax becomes Active
Rule('Bax_to_Bid',
     Bax(b = None, m1 = None,  state='I') + tBid(b = None, s='A')  <>
     Bax(b = 1, m1 = None,  state='I') % tBid(b = 1,  s='A'),
     kbaxbidf, baxbidr)
Rule(Bax(b = 1, m1 = None,  state='I') % tBid(b = 1,  s='A') >>
     Bax(b = None, m1 = None,  state='A') + tBid(b=None, s='A'),
     kbaxbidc)

# Active Bax aggregates into a pore
Rule(Bax(b = None, m = None, state='A') + Bax(b = None, m = None, state='A') <>
     Bax(b = None, m = 1, state='A') % Bax(b = None, m = 1, state='A'),
     kbaxporef, kbaxporer)
Rule(Bax(b = None, m = 1, state='A') % Bax(b = None, m = 1, state='A') >>
     BaxPore(),
     kbaxporec)




     

     
     
