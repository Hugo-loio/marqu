from itertools import product
from functools import reduce

import numpy as np
from scipy.optimize import minimize
import pandas as pd

from .utils import *
from .gauge_transform import *
from .operators import identity_pauli_projector_basis as id_basis
from .operators import pauli_projector_basis

class GeneralizedRateMatrix:
    def __init__(self, nsites : int, local_basis):
        self.nsites = nsites
        basis = [reduce(np.kron, ops).flatten() for ops in 
                 product(local_basis, repeat = nsites)]
        # Change of basis matrix with basis operators as columns
        self.A = np.column_stack(basis)
        self.Adag = np.conj(self.A.T)
        # Pseudo-inverse
        self.Ainv = self.Adag @ np.linalg.inv(self.A @ self.Adag)
        # Lindbladian (matrix form)
        self.lind = np.zeros((4**nsites, 4**nsites), dtype = complex) 
        self.M = np.zeros((len(basis),)*2, dtype = float)
        self.gauge = GeneralizedGaugeTransform(self.A)

    def save(self, name):
        data_dir = check_data_dir()
        rate_matrix_dir = check_dir(data_dir + "rate_matrix/")
        path = check_dir(rate_matrix_dir + name + "/")
        props = {'nsites' : self.nsites}
        np.savetxt(path + "M.csv", self.M + self.gauge.Lambda, 
                   delimiter=',', fmt="%.8f")
        pd.Series(props).to_csv(path + "props.csv", header=False)

    def add_hamiltonian(self, ham):
        eye = np.eye(2**self.nsites)
        self.lind += 1j * (np.kron(ham, eye) - np.kron(eye, ham.T))
        self.M = self._computeM()

    def add_noise(self, jump, rate : float):
        eye = np.eye(2**self.nsites)
        self.lind += rate * np.kron(np.conj(jump.T), jump.T) 
        self.lind -= (rate/2) * np.kron(np.conj(jump.T) @ jump, eye) 
        self.lind -= (rate/2) * np.kron(eye, jump.T @ np.conj(jump)) 
        self.M = self._computeM()

    def _computeM(self):
        res = self.Ainv @ self.lind @ self.A
        if(np.linalg.norm(res.imag) > 1E-10):
            raise RuntimeError("Rate matrix not real")
        return np.real(res).T


class RateMatrix(GeneralizedRateMatrix):
    def __init__(self, nsites : int = 1):
        super().__init__(nsites, pauli_projector_basis(1))

class IdentityLocalRateMatrix(RateMatrix):
    def __init__(self):
        super().__init__(1)
        self.M = np.zeros((7,7))
        self.gauge = IdentityPartialGaugeTransform()

    def add_hamiltonian(self, ham):
        super().add_hamiltonian(ham)
        self.M = np.pad(self.M, ((0,1), (0,1)))
        self.update_identity(mat)

    def add_noise(self, jump, rate : float):
        beforeM = self.M
        super().add_noise(jump, rate)
        beforeM[:-1,:-1] = self.M
        self.M = beforeM

    def update_identity(self, mat):
        for i in range(3):
            val = np.sum(np.abs(mat[:,2*i:2*i+2]))/2
            self.gauge.fixed[-1,2*i:2*i+2] += val
            self.gauge.fixed[2*i:2*i+2,-1] += val/2
            self.gauge.fixed[2*i:2*i+2,2*i:2*i+2] -= val/2
            self.gauge.fixed[-1,-1] -= val
        self.gauge.set(np.zeros(self.gauge.dof))
