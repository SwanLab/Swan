import osqp
import numpy as np
from scipy import sparse
import scipy.sparse as sp

class OsQpSolver:
    def __init__(self, options = None):
        self.options = {'verbose': False,
        'eps_abs': 1e-16, 'eps_rel': 1e-16,
        'eps_prim_inf': 1e-20, 'eps_dual_inf': 1e-20,
                        'polish':1,     
                        'adaptive_rho':False,   #Ensure that result is deterministic
                       # 'scaling':1000,
                        'max_iter':4000}
        if not options is None:
            self.options.update(options)
        self.x = None
        self.y = None
        self.warmstart = False

    def set_tolerance(self, tol_qp):    
        self.options['eps_abs'] = tol_qp
        self.options['eps_rel'] = tol_qp

    def tolerance(self):    
        return self.options['eps_abs']

    def __call__(self,A,dJ,dC,p):   
        n = A.shape[0]
        qtildeEps = dC.shape[0]-p
        zeros = sp.coo_matrix((p+qtildeEps, p+qtildeEps))
        Ps = sp.block_diag((A, zeros))
        qs = np.hstack((dJ, np.zeros((p+qtildeEps,)))).T
        Gs = sp.bmat([[A, -dC.T]])  # AX-DC^T\Lambda = 0
        bound_constraint = sp.bmat(
            [[sp.csc_matrix((qtildeEps, n+p)), sp.eye(qtildeEps, format="csc")]])  # mu>=0
        Gs = sp.vstack((Gs, bound_constraint))
        ls = np.zeros((n+qtildeEps,))
        us = np.zeros((n+qtildeEps,))
        us[n:] = np.inf
        self.warmstart = True
        return self.routine(Ps, qs, Gs, ls, us)

    def routine(self, P, q, G, l=None, u=None):
        P = sp.csc_matrix(P)
        A = sp.csc_matrix(G)
        prob = osqp.OSQP()
        prob.setup(P, q, A, l, u, **self.options)
        if self.warmstart and self.x is not None: 
            try:
                prob.warm_start(x=self.x,y=self.y)
            except ValueError:
                pass
        res = prob.solve()
        self.x = res.x 
        self.y = res.y
        return res.x
