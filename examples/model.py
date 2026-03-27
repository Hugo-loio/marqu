from itertools import product

import numpy as np
import matplotlib.pyplot as plt

from marqu.rate_matrix import RateMatrix
from marqu.gauge_optimizer import lp_optimize

# This is an example script that implements the noisy TFIM used in
# TBA reference to our paper 

paulis = np.array([
    [[0,1],[1,0]],
    [[0,-1j],[1j,0]],
    [[1,0],[0,-1]]
    ])

def plot_mat(mat):
    plt.imshow(mat, aspect='auto') 
    plt.colorbar()
    plt.show()

dim = 2
h = 1
J = 0.5

def local_expand(op, nsites):
    res = np.zeros(np.array(op.shape) * (1 << (nsites-1)), dtype = complex)
    for i in range(nsites):
        res += np.kron(np.eye(1 << i), 
                       np.kron(op, np.eye(1 << (nsites-i-1))))
    return res

# Set the hamiltonian terms
def hamiltonian():
    return model

def generate_uniform(noise):
    # Initialize rate matrix for 2-site interactions
    model  = RateMatrix(2)

    #Set the hamiltonian terms
    model.add_hamiltonian((h/(2*dim))*local_expand(paulis[2], 2))
    model.add_hamiltonian(J*np.kron(paulis[0], paulis[0]))
    
    a = (3/2)*np.abs(J) 
    b = (1/2)*np.abs(J) + np.abs(h)/(2*dim)
    c = (1/2)*np.abs(J) 

    #Set the noise terms
    model.add_noise(np.kron(paulis[0], paulis[0]), a * noise)
    model.add_noise(np.kron(paulis[0], paulis[1]), c * noise)
    model.add_noise(np.kron(paulis[0], paulis[2]), b * noise)
    model.add_noise(np.kron(paulis[1], paulis[0]), c * noise)
    model.add_noise(np.kron(paulis[1], paulis[1]), c * noise)
    model.add_noise(np.kron(paulis[1], paulis[2]), b * noise)
    model.add_noise(np.kron(paulis[2], paulis[0]), b * noise)
    model.add_noise(np.kron(paulis[2], paulis[1]), b * noise)
    model.add_noise(np.kron(paulis[2], paulis[2]), c * noise)

    params, cost = lp_optimize(model)
    model.gauge.set(params)
    plot_mat(model.M + model.gauge.Lambda)
    #model.save(f"tfim_dim{dim:g}_J{J:g}_h{h:g}_noise{noise:g}")

    print("Gamma = ", noise, ", cost:", cost)
    return cost


for noise in np.linspace(0, 1, 4):
    generate_uniform(noise)
