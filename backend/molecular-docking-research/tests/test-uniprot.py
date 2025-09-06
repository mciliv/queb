import unittest

from cure import uniprot


class TestUniprot(unittest.TestCase):
    def test_afabp_get_apo_form(self):
        self.assertEqual("3RZY", uniprot.get_apo_form("AFABP"))

    def test_efabp_get_apo_form(self):
        self.assertEqual("4LKP", uniprot.get_apo_form("EFABP"))
