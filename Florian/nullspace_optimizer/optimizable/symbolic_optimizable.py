from .euclidean_optimizable import EuclideanOptimizable
import sympy as sm
import numpy as np
    
x_sm = sm.IndexedBase('x')

def symbolic_optimizable(optimizable):
    r"""
    A decorator for defining an optimization problem with   
    `Sympy <https://docs.sympy.org/latest/index.html>`_ expressions.    
    Only the objective and constraint functions need to be specified,   
    their derivatives being obtained with   
    symbolic differentiation.   

    **Usage:**
        
    .. code:: python    

       from sympy import sqrt   

       @symbolic_optimizable
       class Problem(EuclideanOptimizable):  
            def x0(self):   
                return [1,1] 
                
            def J(self, x): 
                return x[1]+3/sqrt(x[2])
                
            def H(self, x): 
                return (x[1]*x[0]-2)**2 -1
                    
            # No need to defined dJ and dH which are computed by symbolic   
            # differentiation !
        
    A strong assumption of the decorator is that the optimization variable      
    ``x``
    is a  numpy array or list with a fixed dimension. In the definitions    
    of ``J`` and ``H``, the variable ``x`` is 
    a symbolic    
    `IndexedBase <https://docs.sympy.org/latest/modules/tensor/indexed.html>`_  
    object. The implementation of the class is compatible with the use  
    of symbolic Sympy operators as illustrated above with the square root   
    function ``sqrt``.
    """     
        
    def assign(x):  
        return {x_sm[i]:x[i] for i in range(len(x))}
        
    class SymbolicOptimizable(EuclideanOptimizable): 
        _undecorated = optimizable  
            
        def __init__(self, *args, **kwargs):    
            self._problem = optimizable(*args, **kwargs)
            self.J_sm = self._problem.J(x_sm)
            self.G_sm = self._problem.G(x_sm)
            self.H_sm = self._problem.H(x_sm)
            self.dim = len(self._problem.x0())
            self.diffJ_sm = \
                [self.J_sm.diff(x_sm[i]) for i in range(self.dim)]
            self.diffG_sm = \
                [[Gi.diff(x_sm[i]) for i in range(self.dim)] for Gi in self.G_sm]
            self.diffH_sm = \
                [[Hi.diff(x_sm[i]) for i in range(self.dim)] for Hi in self.H_sm]
            
        def x0(self):   
            return np.asarray(self._problem.x0(),dtype=float)
            
        def J(self, x): 
            return self.J_sm.evalf(subs=assign(x))

        def G(self, x): 
            return np.asarray([Gi.evalf(subs=assign(x)) for Gi in self.G_sm],dtype=float)

        def H(self, x): 
            return np.asarray([Hi.evalf(subs=assign(x)) for Hi in self.H_sm],dtype=float)
            
        def dJ(self, x):    
            return np.asarray([diffJi.evalf(subs=assign(x)) for diffJi in self.diffJ_sm],dtype=float)
            
        def dG(self, x):    
            return np.asarray([ [ diffGi[j].evalf(subs=assign(x)) for j in range(len(diffGi))] for diffGi in self.diffG_sm],dtype=float)

        def dH(self, x):    
            return np.asarray([ [ diffHi[j].evalf(subs=assign(x)) for j in range(len(diffHi))] for diffHi in self.diffH_sm],dtype=float)
                                            
        def accept(self, params, results):  
            self._problem.accept(params, results)

    return SymbolicOptimizable

