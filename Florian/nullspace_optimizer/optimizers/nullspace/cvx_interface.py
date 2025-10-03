import numpy as np
import cvxopt as cvx
from scipy.sparse import find
import scipy.sparse as sp
import scipy.sparse.linalg as lg




class CvxQpSolver:
    """ 
    QP solver using <CVXOPT https://cvxopt.org/index.html>_ algorithm. 
        
    Options can be changed with e.g:

    .. code:: python    

        solver = CvxQpSolver()  
        solver.options['maxiters'] = 200
    """
    def __init__(self, options = None):
        self.options = {'show_progress': False,  
                        'abstol': 1e-20, 
                        'reltol': 1e-20,    
                        'feastol': 1e-20,       
                        'refinement': 3,
                        'maxiters':30}
        self.kktsolver = None
        if not options is None:
            self.options.update(options)
        
    def set_tolerance(self, tol_qp):    
        self.options['abstol'] = tol_qp
        self.options['reltol'] = tol_qp
        self.options['feastol'] = tol_qp

    def tolerance(self):    
        return self.options['abstol']

    def __call__(self,A,dJ,dC,p):   
        n = A.shape[0]
        qtildeEps = dC.shape[0]-p
        zeros = sp.coo_matrix((p+qtildeEps, p+qtildeEps))
        Ps = sp.block_diag((A, zeros))
        qs = np.hstack((dJ, np.zeros((p+qtildeEps,)))).T
        As = sp.bmat([[A, -dC.T]])  # AX-DC^T\Lambda = 0
        bs = np.zeros((n,))
        hs = np.zeros((qtildeEps,))
        Gs = sp.bmat(
            [[sp.csc_matrix((qtildeEps, n+p)),  
              -sp.eye(qtildeEps, format="csc")]])  # mu>=0
        def F(W):
            """A KKT solver for solving 
            [ A    0       A     0  ][ux[:n]]   [bx[:n]]
            [ 0    0     -DC     V  ][ux[n:]]   [bx[n:]]
            [ A  -DC^T     0     0  ][uy]  = [by]
            [ 0    V^T     0   -W^TW][uz]    [bz]
                
            See also <https://cvxopt.org/userguide/coneprog.html#exploiting-structure>_.
            """  
            di = W['di']**2
            S = sp.block_diag((sp.csc_matrix((p,p)),    
                               sp.diags(np.asarray(di).flatten(), format="csc")))
            K = sp.bmat([[A,dC.T],  
                         [dC,-S]],format="csc")
            solve = lg.factorized(K)
            Asolve = lg.factorized(A)
            def f(bx, by, bz):
                rhsK = np.vstack((bx[:n]-by,        
                                  -bx[n:n+p],   
                                  cvx.mul(di,bz)-bx[n+p:]))
                U = solve(rhsK) 
                uy = U[:n]  
                ux_lambda = U[n:]
                ux = np.vstack((Asolve(np.array(bx[:n]))-uy, ux_lambda))
                uz = -np.array(di)*ux_lambda[p:]-cvx.mul(di,bz)
                bx[:] = ux
                by[:] = uy  
                bz[:] = uz*W['d']
                #assert np.allclose(A @ ux[:n]+A@uy,bx[:n])
                #assert np.allclose(-(dC @ uy)[:p],bx[n:][:p])
                #assert np.allclose(-(dC @ uy)[p:]-uz,bx[n:][p:])
                #assert np.allclose(A@ux[:n]-dC.T@ux[n:],by)
                #assert np.allclose(-ux[n:][p:]-W['d']**2*uz,bz)    
                #return (cvx.matrix(ux),cvx.matrix(uy),cvx.matrix(W['d']*uz))
            return f

        I_P, J_P, nz_P = find(Ps)
        I_G, J_G, nz_G = find(Gs)
        I_A, J_A, nz_A = find(As)
        Pcvx = cvx.spmatrix(nz_P, I_P, J_P, Ps.shape) 
        qcvx = cvx.matrix(qs)
        Gcvx = cvx.spmatrix(nz_G, I_G, J_G, Gs.shape)
        hcvx = cvx.matrix(hs)
        Acvx = cvx.spmatrix(nz_A, I_A, J_A, As.shape)
        bcvx = cvx.matrix(bs)
        dims = dict(l=Gs.shape[0],q=[],s=[])
        res = cvx.solvers.coneqp(Pcvx, qcvx, Gcvx, hcvx, dims, Acvx, bcvx,  
                                 initvals = None, kktsolver = F, options=self.options)

        return np.asarray(res['x']).flatten()


    def routine(self, P, q, G, l=None, u=None):   
        print(self.options)
        indicesEq = np.where(l==u)[0]
        A = G[indicesEq,:]  
        b = l[indicesEq]
            
        indicesUp = np.where(np.logical_and(u!=np.inf,l!=u))[0]
        Gup = G[indicesUp,:]
        hup = u[indicesUp]
            
        indicesLow = np.where(np.logical_and(l!=-np.inf,l!=u))[0]    
        Glow = G[indicesLow,:]  
        hlow = l[indicesLow] 
            
        G = sp.vstack((Gup,-Glow))
        h = np.hstack((hup,-hlow))

        I_P, J_P, nz_P = find(P)
        I_G, J_G, nz_G = find(G)
        I_A, J_A, nz_A = find(A)
        Pcvx = cvx.spmatrix(nz_P, I_P, J_P, P.shape) 
        qcvx = cvx.matrix(q)
        Gcvx = cvx.spmatrix(nz_G, I_G, J_G, G.shape)
        hcvx = cvx.matrix(h)
        Acvx = cvx.spmatrix(nz_A, I_A, J_A, A.shape)
        bcvx = cvx.matrix(b)
        dims = dict(l=G.shape[0],q=[],s=[])
            
        kktsolver = self.kktsolver
        res = cvx.solvers.coneqp(Pcvx, qcvx, Gcvx, hcvx, dims, Acvx, bcvx,  
                                 initvals = None, kktsolver = kktsolver, options=self.options)
        return np.asarray(res['x']).flatten()
        
        

