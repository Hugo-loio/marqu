import numpy as np
from scipy.optimize import linprog

from .rate_matrix import RateMatrix

def lp_optimize(rate_matrix : RateMatrix):
    # Constraints Ax <= b, x = [p, t]^T, b = M, A = [-basis -I]
    # Bounds t >= 0
    # Objective function J = sum t 
    mask = ~np.eye(rate_matrix.gauge.dim, dtype = bool) # Off-diagonal mask
    
    b = rate_matrix.M[mask] 
    num_constraints = len(b)
    num_params = rate_matrix.gauge.dof
    
    basis_cols = []
    for k in range(num_params):
        p_basis_vec = np.zeros(num_params)
        p_basis_vec[k] = 1
        rate_matrix.gauge.set(p_basis_vec)
        basis_cols.append(rate_matrix.gauge.Lambda[mask])
        
    basis_mat = np.column_stack(basis_cols)
    A = np.hstack([-basis_mat, -np.eye(num_constraints)])
    
    J = np.concatenate([np.zeros(num_params), np.ones(num_constraints)])
    
    bounds_p = [(None, None) for _ in range(num_params)]
    bounds_t = [(0, None) for _ in range(num_constraints)]
    
    res = linprog(
        c=J,
        A_ub=A,
        b_ub=b,
        bounds=bounds_p + bounds_t,
        method='highs', 
        #options={'disp' : True}
    )
    
    if res.success:
        optimal_p = res.x[:num_params]
        return optimal_p, res.fun
    else:
        raise ValueError(f"Solver failed: {res.message}")
