import numpy as np
import matplotlib.pyplot as plt

from marqu.rate_matrix import PairwiseRateMatrix

ising_M = PairwiseRateMatrix()
ising_M.add_hamiltonian('xx', 1)
#print(np.sum(ising_M.M, axis = 0))
plt.imshow(ising_M.M, aspect='auto')
plt.show()
#print(np.linalg.norm(ising_M.M))
