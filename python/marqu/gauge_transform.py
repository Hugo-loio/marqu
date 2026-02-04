from itertools import product

from .projector import projector

import numpy as np
from scipy.linalg import null_space

def sum_off_diag_negative(matrix):
    negative_matrix = np.where(matrix < 0, matrix, 0)
    np.fill_diagonal(negative_matrix, 0)
    return - np.sum(negative_matrix) 

#class GaugeTransform:
#    def __init__(self, nsites : int):
#        self.nsites = nsites
#        self.dim = 6**nsites
#        self.axis_dim = 3**nsites
#        self.sign_dim = 2**nsites
#        self.line_dof = self.dim - 1
#        self.col_dof = self.axis_dim - 1
#        self.dof = self.line_dof*self.col_dof
#        self.Lambda = np.zeros((self.dim, self.dim))
#
#        axis_indices = np.arange(1, 2*nsites + 1, 2)
#        sign_indices = np.arange(2, 2*nsites + 1, 2)
#        self.perm = np.concatenate(([0], axis_indices, sign_indices))
#        self.inv_perm = np.argsort(self.perm)
#        self.shape_base = ((self.dim, self.dim))
#        self.shape_tensor = (self.dim,) + (3, 2) * self.nsites
#        self.shape_perm_tensor = np.array(self.shape_tensor)[self.perm]
#        self.shape_assign = (self.dim, self.axis_dim, self.sign_dim)
#
#    def set(self, parameters):
#        beta = np.zeros((self.dim, self.axis_dim))
#        beta[:-1,:-1] = np.reshape(parameters, (self.line_dof, self.col_dof))
#        beta[-1,:-1] = -np.sum(beta[:,:-1], axis = 0)
#        beta[:,-1] = -np.sum(beta, axis = 1)
#
#        self.Lambda = np.reshape(self.Lambda, self.shape_tensor)
#        self.Lambda = np.transpose(self.Lambda, self.perm)
#
#        self.Lambda = np.reshape(self.Lambda, self.shape_assign)
#        self.Lambda[:] = beta[:,:,np.newaxis]
#
#        self.Lambda = np.reshape(self.Lambda, self.shape_perm_tensor)
#        self.Lambda = np.transpose(self.Lambda, self.inv_perm)
#
#        self.Lambda = np.reshape(self.Lambda, (self.dim,self.dim))


class GaugeTransform:
    def __init__(self, nsites : int):
        self.nsites = nsites
        self.dim = 6**nsites

        # Configuration to operator space matrix transformation
        A = np.column_stack([projector(i, nsites) for i in range(self.dim)])
        self.basis = null_space(A)
        self.dof = (self.dim-1) * self.basis.shape[1]
        self.Lambda = np.zeros((self.dim, self.dim))

        print("Number of parameters =", self.dof)

    def set(self, parameters):
        for i in range(self.dim-1):
            self.Lambda[i] = self.basis @ parameters[i*self.dim, (i+1)*self.dim]
        self.Lambda[-1] = np.sum(self.Lambda[:-1], axis = 1)
