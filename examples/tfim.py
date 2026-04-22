import sys
from itertools import product

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

from marqu.rate_matrix import RateMatrix
from marqu.gauge_optimizer import lp_optimize
from marqu.operators import  paulis

dim = int(sys.argv[1])
J = float(sys.argv[2])
h = float(sys.argv[3])
gamma = float(sys.argv[4])

def local_expand(op, nsites):
    res = np.zeros(np.array(op.shape) * (1 << (nsites-1)), dtype = complex)
    for i in range(nsites):
        res += np.kron(np.eye(1 << i), 
                       np.kron(op, np.eye(1 << (nsites-i-1))))
    return res

def hamiltonian():
    model  = RateMatrix(2)
    ham = (h/(2*dim))*local_expand(paulis[2], 2) + J*np.kron(paulis[0], paulis[0])  
    model.add_hamiltonian(ham)
    return model

def generate_uniform(gamma):
    model = hamiltonian()
    
    a = (3/2)*np.abs(J) 
    b = (1/2)*np.abs(J) + np.abs(h)/(2*dim)
    c = (1/2)*np.abs(J) 

    model.add_noise(np.kron(paulis[0], paulis[0]), a * gamma)
    model.add_noise(np.kron(paulis[0], paulis[1]), c * gamma)
    model.add_noise(np.kron(paulis[0], paulis[2]), b * gamma)
    model.add_noise(np.kron(paulis[1], paulis[0]), c * gamma)
    model.add_noise(np.kron(paulis[1], paulis[1]), c * gamma)
    model.add_noise(np.kron(paulis[1], paulis[2]), b * gamma)
    model.add_noise(np.kron(paulis[2], paulis[0]), b * gamma)
    model.add_noise(np.kron(paulis[2], paulis[1]), b * gamma)
    model.add_noise(np.kron(paulis[2], paulis[2]), c * gamma)


    params, cost = lp_optimize(model)
    model.gauge.set(params)
    model.save(f"tfim_dim{dim:g}_J{J:g}_h{h:g}_gamma{gamma:g}")

    print("gamma = ", gamma, ", cost = ", cost)

generate_uniform(gamma)

