import unittest

import simplesbml, cobra, pathlib

from OptBoltzmann.maximum_entropy import (
    get_nullspace,
    get_stoichiometric_matrix,
    get_random_initial_fluxes,
    get_standard_change_in_gibbs_free_energy,
    get_initial_beta,
    get_target_log_variable_counts,
    get_fixed_log_counts,
    get_objective_reaction_identifiers,
    get_params,

)

class TestMaximumEntropy(unittest.TestCase):
    """Test that Maximum entropy methods work correctly"""
    def setUp(self):
        """Load test models"""
        self.sbmlfiles =["../models/sbml/split_pathway.xml"]
        self.sbmlmodels = [simplesbml.loadSBMLFile(sbmlfile) for sbmlfile in self.sbmlfiles]
        self.cobramodels = [cobra.io.read_sbml_model(sbmlfile) for sbmlfile in self.sbmlfiles]
    def test_maximum_entropy_pyomo_relaxed(self):
        """Test against analytic solution"""

    def test_get_nullspace(self):
        """Make sure the nullspace is computed correctly"""
        

    def test_get_stoichiometric_matrix(self):
        """Make sure the stoichiometric matrix is extracted from the sbml file correctly"""
        for cobramodel, sbmlmodel in zip(self.cobramodels, self.sbmlmodels):
            expected = cobra.util.array.create_stoichiometric_matrix(cobramodel, array_type='DataFrame')
            actual = get_stoichiometric_matrix(sbmlmodel)
            self.assertTrue(expected.equals(actual))

    def test_get_random_initial_variable_concentrations(self):
        """Make sure the initial variable concentrations pass muster"""

    def test_get_random_initial_fluxes(self):
        """Make sure the fluxes are initialized appropriately"""

    def test_get_standard_change_in_gibbs_free_energy(self):
        """Ensure the gibbs free energy is extracted from the sbml file correctly"""

    def test_get_initial_beta(self):
        """Ensure the initial beta is consistent with the initial fluxes"""

    def test_get_target_log_variable_counts(self):
        """Ensure the log variable counts are extracted from the sbml file correctly"""

    def test_get_fixed_log_counts(self):
        """Ensure the fixed log counts are extracted from the sbml file correctly"""

    def test_get_objective_reaction_identifiers(self):
        """Ensure the flux objective coefficients are extracted from the sbml file correctly"""

    def test_get_params(self):
        """Ensure all parameters are accepted by maximum entropy pyomo"""
