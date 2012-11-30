"""
A suite of tests to ensure that each of the models in the EARM
repository can be successfully loaded and have its reaction network
generated.

For every model, pysb.bng.generate_network(model) is called; if there
are no errors, the test passes.
"""

# Import the full extrinsic apoptosis models:

import earm.lopez_embedded as m1a
import earm.lopez_direct as m2a
import earm.lopez_indirect as m3a
import earm.albeck_11b as m4a
import earm.albeck_11c as m5a
import earm.albeck_11d as m6a
import earm.albeck_11e as m7a
import earm.albeck_11f as m8a
import earm.chen_biophys_j as m9a
import earm.chen_febs_direct as m10a
import earm.chen_febs_indirect as m11a
import earm.cui_direct as m12a
import earm.cui_direct1 as m13a
import earm.cui_direct2 as m14a
import earm.howells as m15a

# Import the MOMP-only models:

import earm.mito.lopez_embedded as m1b
import earm.mito.lopez_direct as m2b
import earm.mito.lopez_indirect as m3b
import earm.mito.albeck_11b as m4b
import earm.mito.albeck_11c as m5b
import earm.mito.albeck_11d as m6b
import earm.mito.albeck_11e as m7b
import earm.mito.albeck_11f as m8b
import earm.mito.chen_biophys_j as m9b
import earm.mito.chen_febs_direct as m10b
import earm.mito.chen_febs_indirect as m11b
import earm.mito.cui_direct as m12b
import earm.mito.cui_direct1 as m13b
import earm.mito.cui_direct2 as m14b
import earm.mito.howells as m15b

from pysb.bng import generate_network
import nose
import traceback

# Tests
# =====

def test_generate_network():
    """Test all models for successful network generation by calling
    :py:func:`check_generate_network` for each model.
    """

    models = [m1a.model, m1b.model,
              m2a.model, m2b.model,
              m3a.model, m3b.model,
              m4a.model, m4b.model,
              m5a.model, m5b.model,
              m6a.model, m6b.model,
              m7a.model, m7b.model,
              m8a.model, m8b.model,
              m9a.model, m9b.model,
              m10a.model, m10b.model,
              m11a.model, m11b.model,
              m12a.model, m12b.model,
              m13a.model, m13b.model,
              m14a.model, m14b.model,
              m15a.model, m15b.model]
    for model in models:
        yield (check_generate_network, model)
    
def check_generate_network(model):
    """Tests that network generation occurs without error for the given
    model."""

    success = False

    try:
        generate_network(model)
        success = True
    except:
        pass

    assert success, "Network generation failed on model %s:\n-----\n%s" % \
                (model.name, traceback.format_exc())

if __name__ == '__main__':
    print "Please run this file as 'nosetests %s'" % __file__

