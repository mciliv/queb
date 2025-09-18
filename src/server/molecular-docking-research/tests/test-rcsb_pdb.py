import unittest

from bio import rcsb_pdb


class TestPDB(unittest.TestCase):
    def test_afabp_apo_id(self):
        self.assertEqual("3RZY", rcsb_pdb.search_protein("FABP4"))

    def test_efabp_apo_id(self):
        self.assertEqual("4LKP", rcsb_pdb.search_protein("FABP5"))


if __name__ == "__main__":
    unittest.main()
