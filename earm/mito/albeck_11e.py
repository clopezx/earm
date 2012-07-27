"""TODO: Docstring"""
from pysb import *
from earm import albeck_modules
from sympy.parsing.sympy_parser import parse_expr

Model()

albeck_modules.momp_monomers()

# The specific MOMP model to use
albeck_modules.albeck_11e()

# Observables
Observable('AcBax_', Bax(bf=None, state='A'))
Observable('cSmac', Smac(state='C'))

