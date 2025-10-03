import numpy as np  
from .euclidean_optimizable import EuclideanOptimizable

def tuple_to_array_optimizable(optimizable):    
    """ 
    A converter for optimization variables containing tuples of arrays to   
    an optimizable with a single array variable.    
        
    This is useful for passing optimization problems with tuples of arrays to IPOPT 
    or MMA solvers.
    """
    class ArrayOptimizable(EuclideanOptimizable):
        _undecorated = optimizable
        def __init__(self, *args, **kwargs): 
            self._problem = optimizable(*args, **kwargs)
            self._problem._results = dict()
            x0 = self._problem.x0()
            n = len(x0)
            indices_float = [True if isinstance(xi, (float,int)) else False for xi in x0]
            x0 = [np.asarray([xi]) if isinstance(xi, (float,int)) else xi for xi in x0]
            cumsums = np.hstack((0,np.cumsum([len(xi) for xi in x0])))
            def array_to_tuple(x):  
                return tuple(x[cumsums[i]] if fl else    
                             x[cumsums[i]:cumsums[i+1]] for (i,fl) in zip(range(n),indices_float))
            self._array_to_tuple = array_to_tuple
            
        def x0(self):   
            return np.hstack(self._problem.x0())
            
        def J(self, x):
            return self._problem.J(self._array_to_tuple(x))

        def G(self, x):
            return self._problem.G(self._array_to_tuple(x))

        def H(self, x):
            return self._problem.H(self._array_to_tuple(x))

        def dJ(self, x):
            return self._problem.dJ(self._array_to_tuple(x))

        def dG(self, x):
            return self._problem.dG(self._array_to_tuple(x))

        def dH(self, x):
            return self._problem.dH(self._array_to_tuple(x))
            
        def retract(self, x, dx):
            newx = self._problem.retract(self._array_to_tuple(x), dx)   
            return np.hstack(newx)
            
        def accept(self, params, results):  
            results_copy = dict()   
            for key in results:
                if key =='x':   
                    results_copy['x'] = [self._array_to_tuple(x) if x is not None   
                                         else None for x in results['x']]  
                else:   
                    results_copy[key] = results[key]
            self._problem.accept(params, results_copy)
     
    return ArrayOptimizable
