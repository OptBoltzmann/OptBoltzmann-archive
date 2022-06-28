import unittest

import simplesbml

from mentos.maximum_entropy import get_params, maximum_entropy_pyomo_relaxed


class TestMaximumEntropy(unittest.TestCase):
    def setUp(self):
        self.model = simplesbml.loadSBMLFile("../models/sbml/split_pathway.xml")

    def test_maximum_entropy_pyomo_relaxed(self):
        maximum_entropy_pyomo_relaxed()

    def test_get_nullspace(self):
        """Make sure the nullspace is computed correctly"""
    def test_get_stoichiometric_matrix(self):
        """Make sure the stoichiometric matrix is extracted from the sbml file correctly"""
    def test_get_random_initial_variable_concentrations(self):
        """Make sure the initial variable concentrations pass muster"""
    def test_get_random_initial_fluxes(self):
        """Make sure the fluxes are initialized appropriately"""
    def test_get_standard_change_in_gibbs_free_energy(self):
        """Ensure the gibbs free energy is extracted from the sbml file correctly"""
    def test_get_initial_beta(self):
        """Ensure the initial beta is consistent with the initial fluxes """
    def test_get_target_log_variable_counts(self):
        """Ensure the log variable counts are extracted from the sbml file correctly"""
    def test_get_fixed_log_counts(self):
        """Ensure the fixed log counts are extracted from the sbml file correctly"""
    def test_get_objective_reaction_identifiers(self):
        """Ensure the flux objective coefficients are extracted from the sbml file correctly"""
    def test_get_params(self):
        """Ensure all parameters are accepted by maximum entropy pyomo"""
