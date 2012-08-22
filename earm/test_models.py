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
import unittest

# Tests
# =====

class TestM1a(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m1a.model)

