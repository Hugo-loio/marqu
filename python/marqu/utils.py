import os

import numpy as np

#Levi-civita tensor
eps = np.zeros((3,3,3), dtype=int)
eps[0,1,2] = 1
eps[1,2,0] = 1
eps[2,0,1] = 1
eps[2,1,0] = -1
eps[0,2,1] = -1
eps[1,0,2] = -1

def pauli_to_int(pauli : str):
    if(pauli == 'x'):
        return 0
    elif(pauli == 'y'):
        return 1
    elif(pauli == 'z'):
        return 2
    else:
        raise ValueError(f"Pauli operator {pauli} not recognized")

def int_to_pauli(alpha : int):
    if(alpha == 0):
        return 'x'
    elif(alpha == 1):
        return 'y'
    elif(alpha == 2):
        return 'z'
    else:
        raise ValueError(f"Pauli operator {alpha} not recognized")


def check_dir(path):
    if(not os.path.isdir(path)):
        try:
            os.mkdir(path)
        except FileExistsError:
            print(path + " was created by another process")
    return path

def check_data_dir():
    path = os.getcwd() + "/marqu_data/"
    return check_dir(path)
