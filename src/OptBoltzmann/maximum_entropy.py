# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 2021

@author: Ethan King
@author: Jeremy Zucker
"""

import itertools

import libsbml
import numpy as np
import numpy.random as nprd
import pandas as pd
import pyomo
import pyomo.environ as pe
import pyutilib.services
import scipy.linalg as spL
import simplesbml
from pyomo.opt import TerminationCondition
from scipy.linalg import norm
from scipy.optimize import least_squares

__all__ = """get_nullspace,
    get_stoichiometric_matrix,
    get_random_initial_fluxes,
    get_standard_change_in_gibbs_free_energy,
    get_initial_beta,
    get_target_log_variable_counts,
    get_fixed_log_counts,
    get_objective_reaction_identifiers,
    get_params""".split(',')

def get_nullspace(A, atol=1e-13, rtol=0):
    """Compute an approximate basis for the nullspace of A.

    The algorithm used by this function is based on the singular value
    decomposition of `A`.

    :param A: ndarray
        A should be at most 2-D.  A 1-D array with length k will be treated
        as a 2-D with shape (1, k)
    :param atol: float
        The absolute tolerance for a zero singular value.  Singular values
        smaller than `atol` are considered to be zero.
    :param rtol : float
        The relative tolerance.  Singular values less than rtol*smax are
        considered to be zero, where smax is the largest singular value.

    If both `atol` and `rtol` are positive, the combined tolerance is the
    maximum of the two; that is::
        tol = max(atol, rtol * smax)
    Singular values smaller than `tol` are considered to be zero.

    Return value
    ------------
    :returns: ns ndarray
        If `A` is an array with shape (m, k), then `ns` will be an array
        with shape (k, n), where n is the estimated dimension of the
        nullspace of `A`.  The columns of `ns` are a basis for the
        nullspace; each element in numpy.dot(A, ns) will be approximately
        zero.
    """

    A = np.atleast_2d(A)
    u, s, vh = svd(A)
    tol = max(atol, rtol * s[0])
    nnz = (s >= tol).sum()
    ns = vh[nnz:].conj().T
    return ns


def get_stoichiometric_matrix(model: simplesbml.SbmlModel):
    """Allocate space for the stoichiometric matrix
    :param model: simplesbml.SbmlModel
    :return stoich: a Dataframe with nrxn columns and nspecies index.
    """
    # Allocate space for the stoichiometry matrix
    stoich = np.zeros((model.getNumSpecies(), model.getNumReactions()))
    for i, speciesId in enumerate(model.getListOfAllSpecies()):
        for j in range(model.getNumReactions()):
            productStoichiometry = 0
            reactantStoichiometry = 0

            numProducts = model.getNumProducts(j)
            for k1 in range(numProducts):
                productId = model.getProduct(j, k1)

                if speciesId == productId:
                    productStoichiometry += model.getProductStoichiometry(j, k1)

            numReactants = model.getNumReactants(j)
            for k1 in range(numReactants):
                reactantId = model.getReactant(j, k1)
                if speciesId == reactantId:
                    reactantStoichiometry += model.getReactantStoichiometry(j, k1)

            st = int(productStoichiometry - reactantStoichiometry)
            stoich[i, j] = st
    return pd.DataFrame(
                stoich, columns=model.getListOfReactions(), index=model.getListOfSpecies()
            )


def get_random_initial_variable_concentrations(model: simplesbml.SbmlModel) -> pd.Series:
    """get vector of initial concnetraions of variable metabolites
    :param model: simplesbml.SbmlModel
    :returns: a random vector of initial concentration of size numFloatingSpecies
    """
    return pd.Series(
        nprd.rand(model.getNumFloatingSpecies()), index=model.getListOfFloatingSpecies()
    )


def get_random_initial_fluxes(model: simplesbml.SbmlModel) -> pd.Series:
    """get vector of initial fluxes
    :param model: simplesbml.SbmlModel
    :returns: a random vector of initial fluxes of size numReactions
    """
    return pd.Series(nprd.rand(model.getNumReactions()), index=model.getListOfReactions())


def get_standard_change_in_gibbs_free_energy(model: simplesbml) -> pd.Series:
    """get standard change in Gibbs free energy
    :param model: simplesbml.SbmlModel
    :returns: pd.Series of change in Gibbs free energy stored as parameters
    """
    return pd.Series(
        [model.getParameterValue(parameter_id) for parameter_id in model.getListOfParameterIds()],
        index=model.getListOfParameterIds(),
    )


def get_initial_beta(nullspace: pd.DataFrame) -> pd.DataFrame:
    r"""get initial beta matrix, where :math:`y=B\beta` for a :math:`\beta\in\mathbb{R}^m`
    :param nullspace: nullspace of stoichiometric matrix
    :returns: random matrix of the size of the nullspace
    """
    return nprd.rand(*nullspace.shape)


def get_target_log_variable_counts(model: simplesbml.SbmlModel) -> pd.Series:
    """get target log variable counts
    :param model: simplesbml.SbmlModel
    :returns: pd.Series target log variable counts extracted from species initial amounts.
    """
    return pd.Series(
        [
            model.getSpeciesInitialAmount(floating_species_id)
            for floating_species_id in model.getListOfFloatingSpecies()
        ],
        index=list(model.getListOfFloatingSpecies()),
    )


def get_fixed_log_counts(model: simplesbml.SbmlModel) -> pd.Series:
    """Get fixed log counts
    :param model: simplesbml.SbmlModel
    :returns: Pandas Series containing inital amounts from boundary species
    """
    return pd.Series(
        [
            model.getSpeciesInitialAmount(fixed_species_id)
            for fixed_species_id in model.getListOfBoundarySpecies()
        ],
        index=list(model.getListOfBoundarySpecies()),
    )


def get_objective_reaction_identifiers(model: simplesbml.SbmlModel) -> pd.Series:
    """Get objective coefficients of reaction identifiers
    :param model: simplesbml.SbmlModel
    :returns: pandas Series containing the reaction ids with nonzero objective coefficients
    """
    model_fbc = model.model.getPlugin("fbc")
    rxn_idx, obj = zip(
        [
            (fluxobj.getReaction(), fluxobj.getCoefficient())
            for fluxobj in model_fbc.getActiveObjective().getListOfFluxObjectives()
        ]
    )
    return pd.Series(obj, index=rxn_idx)


def get_params(model: simplesbml.SbmlModel):
    """
    :param sbml: libsbml.SBMLDocument
    :returns n_ini: vector of initial concentrations of variable metabolites.
    :returns y_ini: initial value of transformed variable used in constraints. length = number of net reactions.
    :returns beta_ini: length = transformed variable used in constraints. length = size of the null space of S.
    :returns target_log_vcounts: vector of experimental (or reasonable estimated) concentrations of variable metabolites.
        These are the constraints for the interior point optimization. (Likewise when we do reinforcement learning,
        these are part of our loss function for learning control)
    :returns f_log_counts: fixed metabolites concentrations
    :returns S: Stoichiometric matrix, n_react x n_M
    :returns K: n_react array of K = exp(-standard change in gibbs free energy/RT)
    :returns obj_rxn_idx: indices of the reactions that are included in the objective function. Currently these are the
        reactions forming proteins, RNA and DNA.
    """
    return dict(
        n_ini=get_random_initial_variable_concentrations(model),
        y_ini=get_random_initial_fluxes(model),
        K=get_standard_change_in_gibbs_free_energy(model),
        S=get_stoichiometric_matrix(model),
        nullspace=get_nullspace(
            S.loc[model.getListOfFloatingSpecies()]
        ),  # Nullspace of nfloatingspecies x nrxns
        beta_ini=get_initial_beta(nullspace),
        target_log_vcounts=get_target_log_variable_counts(model),
        f_log_counts=get_fixed_log_counts(model),
        obj_rxn_idx=get_objective_reaction_identifiers(model),
    )


def maximum_entropy_pyomo_relaxed(
    n_ini, y_ini, beta_ini, target_log_vcounts, f_log_counts, S, K, obj_rxn_idx
):
    """
    :param n_ini: vector of initial concentrations of variable metabolites.
    :param y_ini: initial value of transformed variable used in constraints. length = number of net reactions.
    :param beta_ini: length = transformed variable used in constraints. length = size of the null space of S.
    :param target_log_vcounts: vector of experimental (or reasonable estimated) concentrations of variable metabolites.
        These are the constraints for the interior point optimization. (Likewise when we do reinforcement learning,
        these are part of our loss function for learning control)
    :param f_log_counts: fixed metabolites concentrations
    :param S: Stoichiometric matrix, n_react x n_M
    :param K: n_react array of K = exp(-standard change in gibbs free energy/RT)
    :param obj_rxn_idx: indices of the reactions that are included in the objective function. Currently these are the
        reactions forming proteins, RNA and DNA.
    :returns b_sol: a transformed variable used in constraints = number of reactions
    :returns y_sol: transformed variable used in constraints = number of net reactions
    :returns alpha_sol: reaction regulation
    :returns h_sol: transformed variable used in constraints = number of net reactions
    :returns beta_sol: transformed variable used in constraints = size of null space of S
    :returns n_sol:  vector of solution concentrations of variable metabolites.
    """

    # Flip Stoichiometric Matrix
    S_T = S  # S_T is the Stoich matrix with rows as reactions, columns as metabolites
    S = np.transpose(
        S
    )  # this now is the Stoich matrix with rows metabolites, and columns reactions
    n_react = np.shape(S)[1]

    # Set the System parameters
    FxdM = np.reshape(f_log_counts, (len(f_log_counts), 1))
    K_eq = np.reshape(K, (n_react, 1))

    # Metabolite parms
    n_M_f = len(f_log_counts)  # total number of fixed metabolites
    n_M_v = len(n_ini)  # total number of variable metabolites
    n_M = n_M_f + n_M_v  # total number of metabolites

    ### construct parm indices
    react_idx = np.arange(0, n_react)

    TotM_idx = np.arange(0, n_M)
    VarM_idx = np.arange(0, n_M_v)  # variable metabolite indices
    FxdM_idx = np.arange(n_M_v, n_M)  # fixed metabolite indices

    # Split S into the component corresponding to the variable metabolites S_v
    # and the component corresponding to the fixed metabolites S_f
    S_v_T = np.delete(S_T, FxdM_idx, 1)
    S_v = np.transpose(S_v_T)

    S_f_T = np.delete(S_T, VarM_idx, 1)
    S_f = np.transpose(S_f_T)

    # find a basis for the nullspace of S_v
    Sv_N = spL.null_space(S_v)
    dSv_N = np.shape(Sv_N)[1]  # the dimension of the nullspace

    # Steady State satisfies the following least squares problem
    # || S_T*( SVS* y_hat) - y_hat ||^2 = 0
    # where SVS is the matrix (S_v * S_v_T)^-1 * S_v

    ###precompute SVS
    SvS = np.matmul(np.linalg.inv(np.matmul(S_v, S_v_T)), S_v)

    # Get the system parameters to constuct a model

    # Metabolite parms
    n_M = n_M_f + n_M_v  # total number of metabolites

    ### construct parm indices
    react_idx = np.arange(0, n_react)

    beta_idx = np.arange(0, dSv_N)

    # Set the initial condition
    b_ini = np.matmul(S_v_T, np.reshape(n_ini, (len(n_ini), 1))) + np.matmul(
        S_f_T, np.reshape(FxdM, (n_M_f, 1))
    )
    b_ini = np.reshape(b_ini, len(b_ini))

    h_ini = np.sign(y_ini) * (np.log(2) - np.log(np.abs(y_ini) + np.sqrt(np.power(y_ini, 2) + 4)))

    ## Set the optimization parameters
    VarM_lbnd = -300  # lower bound on the log metabolite counts

    # Pyomo Model Definition
    ######################
    m = pe.ConcreteModel()

    # Input the model parameters

    # set the indices
    #####################
    m.react_idx = pe.Set(initialize=react_idx)

    m.TotM_idx = pe.Set(initialize=TotM_idx)
    m.VarM_idx = pe.Set(initialize=VarM_idx)
    m.FxdM_idx = pe.Set(initialize=FxdM_idx)

    m.beta_idx = pe.Set(initialize=beta_idx)

    m.obj_rxn_idx = pe.Set(initialize=obj_rxn_idx)

    # Stochiometric matrix
    ###########################
    S_idx = list(itertools.product(np.arange(0, n_M), np.arange(0, n_react)))
    S_vals = list(np.reshape(S, [1, n_M * n_react])[0])
    S_dict = dict(list(zip(S_idx, S_vals)))

    m.S = pe.Param(S_idx, initialize=S_dict, mutable=True)

    # variable metabolite stochiometric psuedo inverse matrix
    ############################
    SvS_idx = list(itertools.product(np.arange(0, n_M_v), np.arange(0, n_react)))
    SvS_vals = list(np.reshape(SvS, [1, n_M_v * n_react])[0])
    SvS_dict = dict(list(zip(SvS_idx, SvS_vals)))

    m.SvS = pe.Param(SvS_idx, initialize=SvS_dict)

    ## Nullspace basis Matrix
    ####################
    SvN_idx = list(itertools.product(np.arange(0, n_react), np.arange(0, dSv_N)))
    SvN_vals = list(np.reshape(Sv_N, [1, n_react * dSv_N])[0])
    SvN_dict = dict(list(zip(SvN_idx, SvN_vals)))

    m.SvN = pe.Param(SvN_idx, initialize=SvN_dict)

    # Reaction Equilibrium constants
    ##################
    K_dict = dict(list(zip(react_idx, K)))
    m.K = pe.Param(m.react_idx, initialize=K_dict)

    # Fixed metabolite log counts
    FxdM_dict = dict(list(zip(FxdM_idx, f_log_counts)))
    m.FxdM = pe.Param(m.FxdM_idx, initialize=FxdM_dict)

    # Bounds on the log of the metabolites
    #########
    M_ubnd_dict = dict(list(zip(VarM_idx, target_log_vcounts)))
    m.VarM_ubnd = pe.Param(m.VarM_idx, initialize=M_ubnd_dict)

    # SET the Variables
    #############################

    ## Variable metabolites (log)
    ######################

    Mini_dict = dict(list(zip(VarM_idx, n_ini)))
    m.VarM = pe.Var(VarM_idx, initialize=Mini_dict)

    # steady state fluxes
    yini_dict = dict(list(zip(react_idx, y_ini)))
    m.y = pe.Var(react_idx, initialize=yini_dict)

    # flux null space representation
    betaini_dict = dict(list(zip(beta_idx, beta_ini)))
    m.beta = pe.Var(beta_idx, initialize=betaini_dict)

    # Steady state condition RHS
    bini_dict = dict(list(zip(react_idx, b_ini)))
    m.b = pe.Var(react_idx, initialize=bini_dict)

    hini_dict = dict(list(zip(react_idx, h_ini)))
    m.h = pe.Var(react_idx, initialize=hini_dict)

    # Set the Constraints
    #############################

    def flux_null_rep(m, i):
        """ flux null space representation constraint"""
        return m.y[i] == sum(m.SvN[(i, j)] * m.beta[j] for j in m.beta_idx)

    m.fxnrep_cns = pe.Constraint(m.react_idx, rule=flux_null_rep)

    def steady_state_metab(m, j):
        """steady state metabolism constraint"""
        return m.b[j] == sum(m.S[(k, j)] * m.VarM[k] for k in m.VarM_idx) + sum(
            m.S[(k, j)] * m.FxdM[k] for k in m.FxdM_idx
        )

    m.ssM_cns = pe.Constraint(m.react_idx, rule=steady_state_metab)

    def num_smooth_cns(m, i):
        """Number of smooth constraints"""
        return m.h[i] == (m.y[i] * 1e50 / (abs(m.y[i]) * 1e50 + 1e-50)) * (
            pe.log(2) - pe.log(abs(m.y[i]) + pe.sqrt(m.y[i] ** 2 + 4))
        )

    m.nms_cns = pe.Constraint(m.react_idx, rule=num_smooth_cns)

    # y sign variable
    y_sign_ini_dict = dict(list(zip(react_idx, 0.5 + 0.5 * np.sign(y_ini))))
    m.u = pe.Var(react_idx, bounds=(0, 1), initialize=y_sign_ini_dict)

    Mb = 1000

    def relaxed_reg_cns_upper(m, i):
        """Relaxed upper regulatory constraints"""
        return (m.b[i] - pe.log(m.K[i])) >= m.h[i] - Mb * (m.u[i])

    m.rxr_cns_up = pe.Constraint(m.react_idx, rule=relaxed_reg_cns_upper)

    def relaxed_reg_cns_lower(m, i):
        """Relaxed lower regulatory constraints"""
        return (m.b[i] - pe.log(m.K[i])) <= m.h[i] + Mb * (1 - m.u[i])

    m.rxr_cns_low = pe.Constraint(m.react_idx, rule=relaxed_reg_cns_lower)

    def sign_constraint(m, i):
        """Signed constraints"""
        return (pe.log(m.K[i]) - m.b[i]) * m.y[i] >= 0

    m.sign_y_cns = pe.Constraint(m.react_idx, rule=sign_constraint)

    def y_sign_relax(m, i):
        """Relaxed Flux constraints"""
        return 2 * m.u[i] - 1 == (m.y[i] / (abs(m.y[i]) + 1e-50))

    m.y_sign_relax_cns = pe.Constraint(m.react_idx, rule=y_sign_relax)

    # Variable metabolite upper and lower bounds
    def M_upper_cnstrnts(m, i):
        """Upper constraints on big-M"""
        return m.VarM[i] <= m.VarM_ubnd[i]

    m.VarM_ub_cns = pe.Constraint(m.VarM_idx, rule=M_upper_cnstrnts)

    def M_lower_cnstrnts(m, i):
        """"Lower constraints on big-M"""
        return m.VarM[i] >= VarM_lbnd

    m.VarM_lb_cns = pe.Constraint(m.VarM_idx, rule=M_lower_cnstrnts)

    # Set the Objective function

    # Maximum Entropy Production Objective with subset of reactions obj_rxn_idx

    def _Obj(m):
        """"Set the Objective function

    # Maximum Entropy Production Objective with subset of reactions obj_rxn_idx"""
        # return sum( m.y[j]  *  ( pe.log(m.K[j]) + sum( -m.S[(k,j)]*m.VarM[k] for k in m.VarM_idx ) + sum( -m.S[(w,j)]*m.FxdM[w] for w in m.FxdM_idx )    )    for j in m.obj_rxn_idx )
        ##sum_y = sum( m.y[j]  for j in m.obj_rxn_idx )
        ##return sum( (1 - m.y[j]/sum_y)*m.y[j]  for j in m.obj_rxn_idx )
        return sum(m.y[j] for j in m.obj_rxn_idx)

    m.Obj_fn = pe.Objective(rule=_Obj, sense=pe.maximize)

    # Find a Solution
    ####################

    max_iter = 10000
    max_cpu_time = 800000

    import logging

    logging.getLogger("pyomo.core").setLevel(logging.ERROR)

    # Set the solver to use
    opt = pe.SolverFactory("ipopt", solver_io="nl")

    # Set solver otpions
    opts = {
        "max_iter": max_iter,
        "max_cpu_time": max_cpu_time,
        "tol": 1e-7,
        "acceptable_tol": 1e-6,  # }
        "linear_solver": "ma57",
    }
    #'print_level': 0}
    # 'halt_on_ampl_error': 'yes'}
    # 'dual_inf_tol':1.0,
    # 'acceptable_dual_inf_tol':1.01,
    # 'OF_print_info_string':'yes'}
    #'acceptable_constr_viol_tol':1e-10,
    #'constr_viol_tol': 1e-7,
    #'acceptable_constr_viol_tol':1e-6}
    #'halt_on_ampl_error': 'yes'}

    ## Solve the Model
    status_obj = opt.solve(m, options=opts, tee=True)

    n_sol = np.zeros(n_M_v)
    b_sol = np.zeros(n_react)
    y_sol = np.zeros(n_react)
    beta_sol = np.zeros(dSv_N)
    h_sol = np.zeros(n_react)

    for i in react_idx:
        b_sol[i] = pe.value(m.b[i])
        y_sol[i] = pe.value(m.y[i])
        h_sol[i] = pe.value(m.h[i])

    for i in beta_idx:
        beta_sol[i] = pe.value(m.beta[i])

    for i in VarM_idx:
        n_sol[i] = pe.value(m.VarM[i])

    E_regulation = np.ones(len(y_sol))
    unreg_rxn_flux = np.ravel(rxn_flux(n_sol, f_log_counts, S_T, K, E_regulation))
    alpha_sol = y_sol / unreg_rxn_flux
    alpha_sol = np.ravel(alpha_sol)

    return (b_sol, y_sol, alpha_sol, h_sol, beta_sol, n_sol)


def rxn_flux(v_log_counts, f_log_counts, S, K, E_regulation):
    """Compute initial reaction fluxes"""
    # Flip Stoichiometric Matrix
    S_T = S  # S_T is the Stoich matrix with rows as reactions, columns as metabolites
    S = np.transpose(
        S
    )  # this now is the Stoich matrix with rows metabolites, and columns reactions
    n_react = np.shape(S)[1]
    v_log_counts = np.reshape(v_log_counts, (len(v_log_counts), 1))
    f_log_counts = np.reshape(f_log_counts, (len(f_log_counts), 1))
    tot_log_counts = np.concatenate((v_log_counts, f_log_counts))
    K = np.reshape(K, (len(K), 1))
    E_regulation = np.reshape(E_regulation, (len(E_regulation), 1))

    # forward_odds = K*np.exp(- .25*np.matmul(S_T,tot_log_counts) )*np.exp(- .25*np.matmul(S_T,tot_log_counts) )*np.exp(- .25*np.matmul(S_T,tot_log_counts) )*np.exp(- .25*np.matmul(S_T,tot_log_counts) )
    forward_odds = K * np.exp(-0.25 * np.matmul(S_T, tot_log_counts)) ** 4
    reverse_odds = np.power(K, -1) * np.exp(0.25 * np.matmul(S_T, tot_log_counts)) ** 4

    # forward_odds = np.exp(-.25*(np.matmul(S_T,tot_log_counts) + np.log(K) + np.log(E_regulation) ) )**4
    # reverse_odds = np.exp(.25*(np.matmul(S_T,tot_log_counts) - np.log(K) + np.log(E_regulation) )  )**4

    return E_regulation * (forward_odds - reverse_odds)
