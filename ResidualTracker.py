import numpy as np
class ResidualTracker:
    def __init__(self, A, b):
        self.residuals = []
        self.A = A
        self.b = b

    def __call__(self, xk):
        r = self.b - self.A @ xk
        norm_r = np.linalg.norm(r)
        self.residuals.append(float(norm_r))
