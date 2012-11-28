Tests
=====

The EARM repository contains a number of tests that ensure that the PySB
versions of previously published models accurately duplicate their ODE-based
counterparts. The tests guarantee that this continues to be the case despite
any future changes in PySB or EARM.

Tests are written using the Python testing modules `unittest` and `nose`.

.. toctree::
    :maxdepth: 2

    test_models.rst
    test_albeck_models.rst
    test_shen_models.rst
