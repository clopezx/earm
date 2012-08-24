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

class Test_m1a(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m1a.model)

class Test_m1b(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m1b.model)

class Test_m2a(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m2a.model)

class Test_m2b(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m2b.model)

class Test_m3a(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m3a.model)

class Test_m3b(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m3b.model)

class Test_m4a(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m4a.model)

class Test_m4b(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m4b.model)

class Test_m5a(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m5a.model)

class Test_m5b(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m5b.model)

class Test_m6a(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m6a.model)

class Test_m6b(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m6b.model)

class Test_m7a(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m7a.model)

class Test_m7b(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m7b.model)

class Test_m8a(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m8a.model)

class Test_m8b(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m8b.model)

class Test_m9a(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m9a.model)

class Test_m9b(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m9b.model)

class Test_m10a(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m10a.model)

class Test_m10b(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m10b.model)

class Test_m11a(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m11a.model)

class Test_m11b(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m11b.model)

class Test_m12a(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m12a.model)

class Test_m12b(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m12b.model)

class Test_m13a(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m13a.model)

class Test_m13b(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m13b.model)

class Test_m14a(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m14a.model)

class Test_m14b(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m14b.model)

class Test_m15a(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m15a.model)

class Test_m15b(unittest.TestCase):
    def test_generate_network(self):
        generate_network(m15b.model)

if __name__ == '__main__':
    unittest.main()

