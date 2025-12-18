import numpy as np

def sum_off_diag_negative(matrix):
    negative_matrix = matrix[matrix < 0]
    return - np.sum(matrix[matrix < 0]) + np.sum(np.diagonal(negative_matrix))

class GaugeTransform:
    def __init__(self, nsites : int):
        self.dim1 = 6**nsites - 1
        self.dim2 = 3**nsites - 1
        self.dof = dim1*dim2
        self.Lambda = np.zeros((6**nsites, 6**nsites))

    def set(parameters):
        parameters = np.reshape(parameters, (self.dim1, self.dim2))
        last_row = -np.sum(parameters, axis = 0)
        last_col = -np.sum(parameters, axis = 1)

        self.Lambda[:-1,:-2:2] = parameters
        self.Lambda[:-1,1:-2:2] = parameters
        self.Lambda[-1,:] = last_row
        self.Lambda[:,-2] = last_col
        self.Lambda[:,-1] = last_col
