import numpy as np
import scipy.sparse as sp
from .optimizable import Optimizable

def minmax_optimizable(optimizable):   
    r""" 
    A decorator to implement the min/max formulation for    
    multiobjective optimization.

    It allows to implement the optimization problem     
        
    .. math::
       :label: minmax
    
       \begin{aligned}
           \min_{x\in \mathcal{X}} \max_{1\leqslant k \leqslant r}& \quad J_k(x)\\
           \textrm{s.t.} & \left\{\begin{aligned}
        g_i(x)&=0, \text{ for all } 1\leqslant i\leqslant p,\\
        h_j(x)  &\leqslant  0, \text{ for all }1\leqslant j \leqslant q,\\ 
               \end{aligned}\right.
       \end{aligned}

    by constructing the optimization problem equivalent to      
    :eq:`minmax`: 
    
    .. math::

       \begin{aligned}
           \min_{(x,m)\in \mathcal{X}\times \R} \quad m\\
           \textrm{s.t.} & \left\{\begin{aligned}
        g_i(x)&=0, \text{ for all } 1\leqslant i\leqslant p,\\
        h_j(x)  &\leqslant  0, \text{ for all }1\leqslant j \leqslant q,\\  
        J_k(x) -m & \leqslant 0,\text{ for all }1\leqslant k\leqslant r.
               \end{aligned}\right.
       \end{aligned}
        
    The decorator can be used to implement the multiobjective optimization problem   
    :eq:`minmax` 
    by specifying a multivariate objective function in the definition of the    
    :class:`Optimizable` object as
    follows:

    .. code:: python    
        
       from nullspace_optimizer import EuclideanOptimizable, filtered_optimizable
            
       @minmax_optimizable
       class Problem(EuclideanOptimizable):  
            def x0(self):   
                #x0 
            #...
                
            def J(self, x): 
                # multivariate objective function   
                # return [J1, J2, J3, ..., Jr]
            
            def G(self, x): 
                # equality constraints
                # return [G1,G2,...,Gp] 
        
            def H(self, x): 
                # inequality constraints
                # return [H1,H2,...,Hq] 
    
            def dJ(self, x):    
                # return Jacobian matrix of the multivariate objective  
                # function  
                # return (dJi/dxj)_{i,j}
                
            # etc...
                
    """
        
    class MinMaxOptimizable(Optimizable):
        _undecorated = optimizable

        def __init__(self, *args, **kwargs): 
            self._problem = optimizable(*args, **kwargs) 
            self._problem._results = dict()
            self._Jsave = []
            self._dJsave = []

        def x0(self):   
            old_x0 = self._problem.x0()
            J = self._problem.J(old_x0)
            m0 = max(J)
            return (old_x0,m0)
            
        def J(self, x): 
            return x[1]
            
        def G(self, x): 
            return self._problem.G(x[0])

        def H(self, x): 
            oldH = self._problem.H(x[0])
            J = self._problem.J(x[0])
            return np.hstack((J-x[1],oldH))
            
        def dJ(self, x):    
            old_dJ = sp.csc_matrix(self._problem.dJ(x[0]))
            zeros = np.zeros((1,old_dJ.shape[1]))
            return np.append(zeros,1)
            
        def dG(self, x):    
            old_dG = sp.csc_matrix(self._problem.dG(x[0]))
            if old_dG.shape[0]*old_dG.shape[1] == 0:    
                return []
            zeros = sp.csc_matrix((old_dG.shape[0],1))
            return sp.bmat([[old_dG,zeros]],format="csc")    
            
        def dH(self, x):    
            old_dH = self._problem.dH(x[0])
            zeros = sp.csc_matrix((old_dH.shape[0],1))
                
            old_dJ = self._problem.dJ(x[0])
            minus_one = -np.ones((old_dJ.shape[0],1))
                
            return sp.bmat([[old_dJ,minus_one], 
                            [old_dH,zeros]], format="csc")
            
        def inner_product(self, x): 
            A = self._problem.inner_product(x[0])
            return sp.block_diag((A,1), format="csc")
            
        def retract(self, x, dx):   
            newx = self._problem.retract(x[0], dx[:-1])
            newJ = self._problem.J(newx)
            return (newx, max(newJ))
            
        def accept(self, params, results):    
            (x, m) = results['x'][-1]    
            results_copy = dict()
            J = self._problem.J(x)
            self._Jsave.append(J)

            for key in results:
                if not key in {'x','J','H'}:
                    results_copy[key] = results[key]
                elif key == 'x':    
                    results_copy[key] = [x[0] if x is not None  
                                         else None for x in results['x']]
                elif key == 'J':
                    results_copy[key] = self._Jsave
                elif key == 'H':
                    results_copy[key] = [hi[len(J):] for (hi,J) in zip(results['H'],self._Jsave)]
            params['normalisation_norm'] = lambda xi : np.linalg.norm(xi[:-1],np.inf)
            params['dual_norm'] = lambda dC : sp.csc_matrix.sum(abs(dC[:, :-1]), 1)
                
            self._problem.accept(params, results_copy)
            
    return MinMaxOptimizable
        

