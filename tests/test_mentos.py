import unittest
import libsbml
from  mentos.maximum_entropy import maximum_entropy_pyomo_relaxed, get_params
class TestMaximumEntropy(unittest.TestCase):

    def setUp(self):
        r = libsbml.SBMReader()
        sbml = r.readSBML("../models/sbml/split_pathway.xml")
        n_ini, y_ini, beta_ini, target_log_vcounts, f_log_counts, S, K, obj_rxn_idx = get_params(sbml)

    def test_maximum_entropy_pyomo_relaxed(self):
        maximum_entropy_pyomo_relaxed()