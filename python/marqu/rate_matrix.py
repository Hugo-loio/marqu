from abc import ABC, abstractmethod

import numpy as np
from scipy.optimize import minimize
import pandas as pd

from .utils import *
from .gauge_transform import *
from .rate_matrix_terms import *

#def cost_function(params, rate_matrix, alpha):
#    rate_matrix.gauge.set(params)
#    matrix = rate_matrix.M + rate_matrix.gauge.Lambda
#    #negative_matrix = np.where(matrix < 0, matrix, 0)
#    np.fill_diagonal(matrix, 0)
#
#    # Check for overflow: if -x * alpha is too large, exp() overflows.
#    # For very large negative numbers, softplus(-x) ~= -x.
#    # We can use np.logaddexp for numerical stability: log(exp(a) + exp(b))
#    
#    # We want sum( log(1 + exp(-alpha * x)) ) / alpha
#    # This is equivalent to sum( logaddexp(0, -alpha * x) ) / alpha
#    return np.sum(np.logaddexp(0, -alpha * matrix)) / alpha

# Add the noise 
class RateMatrix:
    def __init__(self, nsites : int = 1):
        self.nsites = nsites
        self.M = np.zeros((6**nsites,6**nsites))
        self.gauge = GaugeTransform(nsites)

    def save(self, name):
        data_dir = check_data_dir()
        rate_matrix_dir = check_dir(data_dir + "rate_matrix/")
        path = check_dir(rate_matrix_dir + name + "/")
        props = {'nsites' : self.nsites}
        np.savetxt(path + "M.csv", self.M, delimiter=',', fmt="%.8f")
        pd.Series(props).to_csv(path + "props.csv", header=False)

    @abstractmethod
    def add_hamiltonian(self, pauli : str, coupling : complex):
        raise NotImplementedError

    #def gauge_optimize(self):
    #    zeros = np.zeros(self.gauge.dof)
    #    params = np.copy(zeros)
    #    print("No gauge cost: ", sum_off_diag_negative(self.M))
    #    alphas = [100.0, 1000.0, 10000.0]
    #    for alpha in alphas:
    #        result = minimize(
    #                        fun=cost_function, 
    #                        x0=params, 
    #                        args=(self, alpha), 
    #                        method='L-BFGS-B',
    #                        options={
    #                            'disp': True, 
    #                            'maxiter': 5000, 
    #                            'maxfun': 100000,
    #                            'ftol': 1e-12, 
    #                            'gtol': 1e-12
    #                        }
    #                    )

    #        print(result)
    #        self.gauge.set(result.x)
    #        params = result.x
    #        print("New cost: ", sum_off_diag_negative(self.M + self.gauge.Lambda))

class LocalRateMatrix(RateMatrix):
    def __init__(self):
        super().__init__(1)

    def add_hamiltonian(self, pauli : str, coupling : complex):
        self.M += local_hamiltonian_M(pauli_to_int(pauli), coupling)

    def add_noise(self, a : NDArray[complex], b : complex, damping_rate : float):
        self.M += damping_rate*local_noise_M(a,b)

class PairwiseRateMatrix(RateMatrix):
    def __init__(self):
        super().__init__(2)

    # Identity is a blank space in the pauli, ex: pauli = ' x'
    def add_hamiltonian(self, pauli : str, coupling : complex):
        if(pauli[0] == ' '):
            pauli2 = pauli_to_int(pauli[1])
            self.M += np.kron(np.eye(6), local_hamiltonian_M(pauli2, coupling))
        elif(pauli[1] == ' '):
            pauli1 = pauli_to_int(pauli[0])
            self.M += np.kron(local_hamiltonian_M(pauli1, coupling), np.eye(6))
        else:
            pauli1 = pauli_to_int(pauli[0])
            pauli2 = pauli_to_int(pauli[1])
            self.M += pairwise_hamiltonian_M(pauli1, pauli2, coupling)

    def add_local_noise(self, a : NDArray[complex], b : complex, damping_rate : float):
        auxM = damping_rate*local_noise_M(a,b)
        self.M += np.kron(np.eye(6), auxM)
        self.M += np.kron(auxM, np.eye(6))
    
    def add_pauli_pair_noise(self, pauli : str, dampling_rate : float):
        pauli1 = pauli_to_int(pauli[0])
        pauli2 = pauli_to_int(pauli[1])
        self.M += damping_rate * pauli_pair_noise_M(pauli1, pauli2)
