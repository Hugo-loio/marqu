from itertools import product

import numpy as np
from numpy.typing import NDArray

from .utils import *

def flatten(axis, sign):
    return axis*2 + (1 - sign) // 2

def g_tensor(aj, ak, pj1, sj1, pk1, sk1, pj2, sj2, pk2, sk2):
    res1 = eps[aj,pj1,pj2]*sj1*sj2*(sk2 + sk1*int(ak == pk1))*int(ak == pk2)
    res2 = eps[ak,pk1,pk2]*sk1*sk2*(sj2 + sj1*int(aj == pj1))*int(aj == pj2)
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
    a = np.array(a, dtype = complex)

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

def pauli_pair_noise_M(pauli1 : int, pauli2 : int):
    paulis = [0,1,2]
    signs = [-1,1]
    M = np.zeros((36,36))
    for pj,sj,pk,sk in product(*((paulis, signs,)*2)):
        row = 6*flatten(pj, sj) + flatten(pk, sk)
        M[row, row] -= 1 
        col = 6*flatten(pj, (2*int(pauli1 == pj) - 1)*sj) + \
                flatten(pk, (2*int(pauli2 == pk) - 1)*sk) 
        M[row, col] += 1
    return M
