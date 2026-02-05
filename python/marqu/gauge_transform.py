from itertools import product

from .projector import projector

import numpy as np
from scipy.linalg import null_space

def sum_off_diag_negative(matrix):
    negative_matrix = np.where(matrix < 0, matrix, 0)
    np.fill_diagonal(negative_matrix, 0)
    return - np.sum(negative_matrix) 

class GaugeTransform:
    def __init__(self, nsites : int):
        self.nsites = nsites
        self.dim = 6**nsites

        # Configuration to operator space matrix transformation
        A = np.column_stack([projector(i, nsites).flatten()
                             for i in range(self.dim)])
        A = np.vstack([A.real, A.imag])
        #print(A.shape)
        self.basis = null_space(A)
        self.dof = (self.dim-1) * self.basis.shape[1]
        self.row_dof = self.basis.shape[1]
        self.Lambda = np.zeros((self.dim, self.dim))
        #print(self.basis)
        #print(self.basis.shape, self.Lambda.shape, self.dof, self.dim)

        #print("Number of parameters =", self.dof)

    def set(self, parameters):
        for i in range(self.dim-1):
            self.Lambda[i] = self.basis @ \
                    parameters[i*self.row_dof:(i+1)*self.row_dof]
        self.Lambda[-1] = np.sum(self.Lambda[:-1], axis = 0)
