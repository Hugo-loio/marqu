import numpy as np

paulis = np.array([
    [[0,1],[1,0]],
    [[0,-1j],[1j,0]],
    [[1,0],[0,-1]]
    ])


# Many-body projector operator for configuration index C 
def projector(C : int, nsites : int):
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
