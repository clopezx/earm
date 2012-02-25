from pysb import *

Model()

Parameter('A', 1)
Parameter('I', 1)

Monomer('X', ['c', 'd'])
Monomer('Y', 'e')

Rule('xxy', X(c=1, d=None) % X(c=1, d=None) + Y(e=None) >> X(c=1, d=2) % X(c=1, d=3) % Y(e=[2,3]), A)

Initial(X(c=1, d=None) % X(c=1, d=None), I)
Initial(Y(e=None), I)
