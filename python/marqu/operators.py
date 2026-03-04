from itertools import product
from functools import reduce

import numpy as np

paulis = np.array([
    [[0,1],[1,0]],
    [[0,-1j],[1j,0]],
    [[1,0],[0,-1]]
    ])


# Many-body projector operator for configuration index C 
def pauli_projector(C : int, nsites : int):
    axes = np.zeros(nsites, dtype = int)
    signs = np.zeros(nsites, dtype = int)
    for i in range(nsites - 1, -1, -1):
        axes[i] = (C % 6) // 2
        signs[i] = 1 - 2 * (C % 2)
        C //= 6

    projector = 1 
    for i in range(nsites):
        projector = np.kron(projector, (np.eye(2) + signs[i]*paulis[axes[i]])/2)
    return projector


def pauli_projector_basis(nsites : int):
    return [pauli_projector(i, nsites) for i in range(6**nsites)]


def identity_pauli_projector_basis(nsites : int):
    local_basis = [pauli_projector(i, 1) for i in range(6)] + [np.eye(2)]
    return [reduce(np.kron, ops) for ops in product(local_basis, repeat=nsites)]
