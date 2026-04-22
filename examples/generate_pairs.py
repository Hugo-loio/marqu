import sys
import os
from itertools import product
import random

import numpy as np

random.seed(1)

N = int(sys.argv[1])
dim = int(sys.argv[2])

sites = list(product(np.arange(N), repeat = dim))

shift = [N//2] * dim
center_site = tuple(np.array(sites[0]) + shift)
pairs = [[sites[1], sites[2]], [sites[0], center_site]]
sites.remove(center_site)
sites.remove(sites[2])
sites.remove(sites[1])
sites.remove(sites[0])

for _ in range(len(sites)//2):
    site1 = random.choice(sites)
    sites.remove(site1)
    site2 = random.choice(sites)
    sites.remove(site2)
    pairs.append([site1, site2])

pairs = np.array(pairs)

flatten_rescale = np.array([N**i for i in range(dim)][::-1])
flat_pairs = pairs * flatten_rescale[np.newaxis, np.newaxis, :]
flat_pairs = np.sum(flat_pairs, axis = 2)

print(flat_pairs)

os.makedirs("pairs", exist_ok=True)
np.savetxt(f"pairs/dim{dim}_N{N}.csv", flat_pairs, fmt='%d', delimiter=',')
