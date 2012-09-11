"""
EARM: Extrinsic Apoptosis Reaction Model
========================================

EARM is a family of models of the extrinsic apoptosis pathway written using the
Python modeling framework PySB. Model variants focus on exploring alternative
hypotheses for the mechanism of Mitochondrial Outer Membrane Permeabilization
(MOMP), controlled by the Bcl-2 family of proteins.

Documentation is available in the docstrings, doc/ directory and online at
http://earm.readthedocs.org .

Modules
-------
::

 albeck_modules   --- components for albeck_* models
 lopez_modules    --- components for lopez_* models
 shen_modules     --- components for chen_*, cui_* and howells models
 shared           --- shared constants and macros for all models

 everything else (including mito.*)
                  --- the models

Scripts
-------
::

 estimate.py      --- simple parameter estimation using simulated annealing
 model_specs.py   --- display number of rules/odes/params for all models
 test_models.py   --- minimal test to ensure models contain no blatant errors

"""
