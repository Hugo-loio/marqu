from itertools import product
from abc import ABC, abstractmethod

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from .utils import *
from .gauge_transform import *

def flatten(axis, sign):
    return axis*2 + (1 - sign) // 2

def g_tensor(aj, ak, pj1, sj1, pk1, sk1, pj2, sj2, pk2, sk2):
    res1 = eps[aj,pj1,pj2]*sj1*sj2*(sk2 + sk1*int(ak == pk1))*int(ak == pk2)
    res2 = eps[ak,pk1,pk2]*sk1*sk2*(sj2 + sj1*int(aj == pj1))*int(aj == pj2)
    print(eps[aj,pj1,pj2], sj1, sj2, ak, pk2)
    return res1 + res2

def local_hamiltonian_M(pauli : int, coupling : float):
    paulis = [0,1,2]
    signs = [-1,1]
    M = np.zeros((6,6))
    for pi,si,pj,sj in product(paulis, signs, paulis, signs):
        i = pi*2 + (si + 1)//2
        j = pj*2 + (sj + 1)//2
        M[i,j] = -coupling*eps[pauli, pi, pj]*si*sj
    return M

def pairwise_hamiltonian_M(pauli1 : int, pauli2 : int, coupling : complex):
    paulis = [0,1,2]
    signs = [-1,1]
    M = np.zeros((36,36))
    for pj1,sj1,pk1,sk1,pj2,sj2,pk2,sk2 in product(*((paulis, signs,)*4)):
        row = 6*flatten(pj1, sj1) + flatten(pk1, sk1)
        col = 6*flatten(pj2, sj2) + flatten(pk2, sk2)
        M[row, col] = -0.5*coupling*g_tensor(pauli1, pauli2, 
                                             pj1, sj1, pk1, sk1, 
                                             pj2, sj2, pk2, sk2)
    return M

def local_traceless_noise_M(a : NDArray[complex]):
    M = np.zeros((6,6))
    a2 = np.real(a * np.conjugate(a))

    paulis = [0,1,2]
    signs = [-1,1]
    #gamma = beta_j', s = s_j'
    for betaj,sj,gamma,s in product(*((paulis, signs,)*2)):
        row = flatten(betaj, sj)
        col = flatten(gamma, s)
        aux1 = np.sum([np.imag(a[alpha] * np.conjugate(a[gamma])) * eps[alpha, gamma, betaj] for alpha in paulis for gamma in paulis])
        aux2 = a2[betaj] - np.sum(a2)
        if(gamma == betaj): 
            M[row, col] += sj * aux1
            M[row, col] += sj * s * aux2
        else:
            M[row, col] += sj * s * np.real(a[betaj] * np.conjugate(a[gamma]))

    return M

def local_noise_M(a : NDArray[complex], b : complex):
    M = local_traceless_noise_M(a)

    paulis = [0,1,2]
    signs = [-1,1]
    #gamma = beta_j', s = s_j'
    for betaj,sj,gamma,s in product(*((paulis, signs,)*2)):
        row = flatten(betaj, sj)
        col = flatten(gamma, s)
        for alpha in paulis:
            M[row, col] += np.imag(b*np.conjugate(a[alpha])) * \
                    eps[alpha, betaj, gamma] * sj * s

    return M

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

    #@abstractmethod
    #def add_noise(self, pauli : str, dampling_rate : float):
    #    raise NotImplementedError

    @abstractmethod
    def add_hamiltonian(self, pauli : str, coupling : complex):
        raise NotImplementedError

    def gauge_optimize(self):
        params = np.zeros(self.gauge.dof)



class LocalRateMatrix(RateMatrix):
    def __init__(self, name):
        super().__init__(1)

    def add_hamiltonian(self, pauli : str, coupling : complex):
        self.M += local_hamiltonian_M(pauli_to_int(pauli), coupling)


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
