from itertools import product

from .operators import pauli_projector_basis

import numpy as np
from scipy.linalg import null_space

def sum_off_diag_negative(matrix):
    negative_matrix = np.where(matrix < 0, matrix, 0)
    np.fill_diagonal(negative_matrix, 0)
    return - np.sum(negative_matrix) 

class GeneralizedGaugeTransform:
    def __init__(self, nsites : int, operator_basis):
        self.nsites = nsites

        op_basis = operator_basis(nsites)
        A = np.column_stack([op.flatten() for op in op_basis])
        A = np.vstack([A.real, A.imag])

        self.dim = len(op_basis)
        self.basis = null_space(A)
        self.dof = (self.dim-1) * self.basis.shape[1]
        self.row_dof = self.basis.shape[1]
        self.Lambda = np.zeros((self.dim, self.dim))

    def set(self, parameters):
        for i in range(self.dim-1):
            self.Lambda[i] = self.basis @ \
                    parameters[i*self.row_dof:(i+1)*self.row_dof]
        self.Lambda[-1] = -np.sum(self.Lambda[:-1], axis = 0)

class GaugeTransform(GeneralizedGaugeTransform):
    def __init__(self, nsites : int):
        super().__init__(nsites, pauli_projector_basis)

# One site, extend to the identity configuration without applying a gauge
class IdentityPartialGaugeTransform(GaugeTransform):
    def __init__(self):
        super().__init__(1)
        self.dim += 1
        self.Lambda = np.pad(self.Lambda, ((0,1), (0,1)))
        self.fixed = np.zeros(self.Lambda.shape)

    def set(self, parameters):
        self.Lambda[:] = 0
        for i in range(self.dim-2):
            self.Lambda[i,:-1] = self.basis @ \
                    parameters[i*self.row_dof:(i+1)*self.row_dof]
        self.Lambda[-2,:-1] = - np.sum(self.Lambda[:-2,:-1], axis=0) 
        self.Lambda += self.fixed
