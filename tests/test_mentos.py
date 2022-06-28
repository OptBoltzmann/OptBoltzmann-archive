import unittest

import simplesbml

from mentos.maximum_entropy import get_params, maximum_entropy_pyomo_relaxed


class TestMaximumEntropy(unittest.TestCase):
    def setUp(self):
        self.model = simplesbml.loadSBMLFile("../models/sbml/split_pathway.xml")

    def test_maximum_entropy_pyomo_relaxed(self):
        maximum_entropy_pyomo_relaxed()

    def test_get_nullspace(self):
        