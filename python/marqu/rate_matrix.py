from itertools import product
from abc import ABC, abstractmethod

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from .utils import eps
from .utils import pauli_to_int
from .utils import check_data_dir
from .utils import check_dir

def g_tensor(aj, ak, pj1, sj1, pk1, sk1, pj2, sj2, pk2, sk2):
    res1 = eps[aj,pj1,pj2]*sj1*sj2*(sk2 + sk1*int(ak == pk1))*int(ak == pk2)
    res2 = eps[ak,pk1,pk2]*sk1*sk2*(sj2 + sj1*int(aj == pj1))*int(aj == pj2)
    print(eps[aj,pj1,pj2], sj1, sj2, ak, pk2)
    return res1 + res2

def all_pairs(N : int, pbc : bool, distance : int = 1):
    return [(i, (i + distance) % N) for i in range(N) 
            if (i + distance < N or pbc)]

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
        row = 6*(pj1*2 + (sj1 + 1)//2) + pk1*2 + (sk1 + 1)//2
        col = 6*(pj2*2 + (sj2 + 1)//2) + pk2*2 + (sk2 + 1)//2
        M[row, col] = -0.5*coupling*g_tensor(pauli1, pauli2, 
                                             pj1, sj1, pk1, sk1, 
                                             pj2, sj2, pk2, sk2)
        if(M[row, col] != 0):
            print(M[row, col])
    return M

def local_traceless_noise_M(a : NDArray[complex]):
    paulis = [0,1,2]
    signs = [-1,1]
    M = np.zeros((6,6))
    a2 = np.square(a)
    px = np.array([[0,1],[1,0]])

    M += kron(np.diag(a2 - 1), np.eye(2))
    M += kron(np.diag(np.sum(a2) - np.square(a)), px)
    M += kron(np.diag(np.real(np.conjugate(a)*np.sum(a) - a2)), np.eye(2)-px)

    return M


# Add the noise 
class RateMatrix:
    def __init__(self, nsites : int = 1):
        self.nsites = nsites
        self.M = np.zeros((6**nsites,6**nsites))

    @abstractmethod
    def add_hamiltonian(self, pauli : str, coupling : complex):
        raise NotImplementedError

    #@abstractmethod
    #def add_noise(self, pauli : str, dampling_rate : float):
    #    raise NotImplementedError

    def save(self, name):
        data_dir = check_data_dir()
        rate_matrix_dir = check_dir(data_dir + "rate_matrix/")
        path = check_dir(rate_matrix_dir + name + "/")
        props = {'nsites' : self.nsites}
        np.savetxt(path + "M.csv", self.M, delimiter=',', fmt="%.8f")
        pd.Series(props).to_csv(path + "props.csv", header=False)

    #def gauge_optimize(self):

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
