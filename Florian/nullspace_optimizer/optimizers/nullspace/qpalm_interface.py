import qpalm
import numpy as np
import scipy.sparse as sp
 

class QPALM_Solver:
    
    def __init__(self, options = None):
        self.options = {'verbose': False,
        'eps_abs': 1e-14, 'eps_rel': 1e-14,
        'eps_prim_inf': 1e-20, 'eps_dual_inf': 1e-20,   
        'max_iter': 400
        }
        if not options is None:
            self.options.update(options)
        self.x = None
        self_Y = None
        self.save = True
    
    def set_tolerance(self, tol_qp):
        self.options['eps_abs'] = tol_qp
    
    def tolerance(self):
        return self.options['eps_abs']
    
    def __call__(self, A, dJ, dC, p):
        n = A.shape[0]
        qtildeEps = dC.shape[0]-p
        zeros = sp.coo_matrix((p+qtildeEps, p+qtildeEps))
        Ps = sp.block_diag((A, zeros))
        qs = np.hstack((dJ, np.zeros((p+qtildeEps,)))).T
        Gs = sp.bmat([[A, -dC.T]])  # AX-DC^T\Lambda = 0
        bound_constraint = sp.bmat([[sp.csc_matrix((qtildeEps, n+p)), sp.eye(qtildeEps, format="csc")]])  # mu>=0
        Gs = sp.vstack((Gs, bound_constraint))
        ls = np.zeros((n+qtildeEps,))
        us = np.zeros((n+qtildeEps,))
        us[n:] = np.inf
        return self.routine(Ps, qs, Gs, ls, us) 

    def routine(self, P, q, G, l=None, u=None):
        
        # Size of the problem
        data = qpalm.Data(P.shape[0], G.shape[0])
        
        data.Q = sp.csc_matrix(P) 
        data.A = sp.csc_matrix(G)
        data.q = q
        data.bmin = l
        data.bmax = u
    
        # %% Configure the solver
        settings = qpalm.Settings()
        for setting, val in self.options.items():   
            exec(f"settings.{setting} = val")
            
    
        # %% Solve the problem
        solver = qpalm.Solver(data, settings)
        solver.solve()
    
        # %% Print the results
        #print("Status:     ", solver.info.status)
        #print("Solution:   ", solver.solution.x)
        #print("Multipliers:", solver.solution.y)
      
        self.x = solver.solution.x
        self.y = solver.solution.y
        
        return solver.solution.x
