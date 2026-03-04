from abc import ABC, abstractmethod

import numpy as np
from scipy.optimize import minimize
import pandas as pd

from .utils import *
from .gauge_transform import *
from .rate_matrix_terms import *
from .operators import identity_pauli_projector_basis as id_basis

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
        np.savetxt(path + "M.csv", self.M + self.gauge.Lambda, 
                   delimiter=',', fmt="%.8f")
        pd.Series(props).to_csv(path + "props.csv", header=False)

    @abstractmethod
    def add_hamiltonian(self, pauli : str, coupling : complex):
        raise NotImplementedError


class LocalRateMatrix(RateMatrix):
    def __init__(self):
        super().__init__(1)

    def add_hamiltonian(self, pauli : str, coupling : complex):
        self.M += local_hamiltonian_M(pauli_to_int(pauli), coupling)

    def add_noise(self, a : NDArray[complex], b : complex, damping_rate : float):
        self.M += damping_rate*local_noise_M(a,b)
        

class IdentityLocalRateMatrix(RateMatrix):
    def __init__(self):
        self.nsites = 1
        self.M = np.zeros((7,7))
        self.gauge = IdentityPartialGaugeTransform()

    def add_hamiltonian(self, pauli : str, coupling : complex):
        mat = local_hamiltonian_M(pauli_to_int(pauli), coupling)
        self.M[:-1,:-1] += mat
        self.update_identity(mat)

    def add_noise(self, a : NDArray[complex], b : complex, damping_rate : float):
        self.M[:-1,:-1] += damping_rate*local_noise_M(a,b)

    def update_identity(self, mat):
        for i in range(3):
            val = np.sum(np.abs(mat[:,2*i:2*i+2]))/2
            self.gauge.fixed[-1,2*i:2*i+2] += val
            self.gauge.fixed[2*i:2*i+2,-1] += val/2
            self.gauge.fixed[2*i:2*i+2,2*i:2*i+2] -= val/2
            self.gauge.fixed[-1,-1] -= val
            #print(self.gauge.fixed)
        #self.gauge.fixed[:] = 0
        self.gauge.set(np.zeros(self.gauge.dof))


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
    
    def add_pauli_pair_noise(self, pauli : str, damping_rate : float):
        pauli1 = pauli_to_int(pauli[0])
        pauli2 = pauli_to_int(pauli[1])
        self.M += damping_rate * pauli_pair_noise_M(pauli1, pauli2)


'''
class IdentityPairwiseRateMatrix(RateMatrix):
    def __init__(self):
        self.nsites = 2
        self.gauge = GeneralizedGaugeTransform(2, id_basis)
        self.M = np.zeros((self.gauge.dim,)*2)

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
    
    def add_pauli_pair_noise(self, pauli : str, damping_rate : float):
        pauli1 = pauli_to_int(pauli[0])
        pauli2 = pauli_to_int(pauli[1])
        self.M += damping_rate * pauli_pair_noise_M(pauli1, pauli2)
        '''
