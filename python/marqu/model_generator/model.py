import numpy as np


class Model:
    def __init__(self, name : str):
        self.name = name
        self.is_ti = True
        self.ham_pairwise = []
        self.ham_local = []
        self.noise_local = []
