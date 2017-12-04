import unittest
import sys
import load_data
sys.path.append('commons')
from common_utils import *

class utilsTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_clean_variant(self):
        variants = {
                '4-15987755-G-C':'4-15987755-G-C',
                '4-15987755-GC-CC':'4-15987755-G-C',
                '4-15987755---C':'4-15987754-T-TC',
                }
        for k,v in variants.items():
            self.assertEqual(clean_variant(k),v)
        

if __name__ == '__main__':
    unittest.main()
