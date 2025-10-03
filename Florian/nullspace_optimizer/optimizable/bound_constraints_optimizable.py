import numpy as np
import scipy.sparse as sp
from .optimizable import Optimizable

def bound_constraints_optimizable(l=None,u=None): 
    r"""
    A decorator for adding bound constraints to an optimizable object.
        
    **Usage:**
        
    .. code:: python    
            
       @bound_constraints_optimizable(l,u)
       class Problem(Optimizable):  
            def x0(self):   
                #x0 
            #...

    This decorator appends bound constraints ``l<=x<=u``    
    to an :class:`Optimizable` object of the form  

    .. math::
    
       \begin{aligned}
           \min_{x\in \R^n}&  \quad J(x)\\
           \textrm{s.t.} & \left\{\begin{aligned}
        g(x)&=0,\\
        h(x)&\leqslant 0,
               \end{aligned}\right.
       \end{aligned}

    The optimization variable ``x`` can be: 
        
    * a numpy array, in that case   
      ``l`` and ``u`` need to be either:    

      * a floating number (will implement ``l<=x<=u``)
      * a numpy arrays of the same size of ``x``. 
        (will implement ``l[i]<=x[i]<=u[i]`` for all ``i``)
      * ``None`` in order to prescribe only a lower or an upper bound    
        (for instance, if ``l`` is ``None``, then the decorator implements   
        ``x<=u``).

    * a tuple of numpy arrays:
      ``x = (x_1, x_2, ..., x_n)`` with fixed dimension ``n``.  
      In that case, ``l`` and ``u`` need to tuples of lower and upper bounds:

      .. code:: python

         l = (l_1, l_2, ..., l_n) # lower bounds 
         u = (u_1, u_2, ..., u_n) # upper bounds

      The decorator implements then the bound constraints ``l_i<=x_i<=u_i``   
      for every ``1<=i<=n`` with the rule above if ``l_i`` and ``u_i`` are    
      floating numbers or numpy arrays. 
       
    .. note::

       If both ``l_i`` and ``u_i`` are ``None``, then no restriction  
       applies to the optimization variable ``x_i``  which
       needs not be a numpy array.
    """
    if not isinstance(l, tuple):
        l = (l,)
    else:
        l = l 
    if not isinstance(u, tuple):
        u = (u,)
    else:
        u = u

    def bound_constraints_optimizable_impl(optimizable):
        class BoundConstraintOptimizable(Optimizable):
            _undecorated = optimizable

            def __init__(self, *args, **kwargs):    
                self._problem = optimizable(*args, **kwargs)

            def x0(self):   
                return self._problem.x0()

            def J(self, x): 
                return self._problem.J(x)   

            def dJ(self, x): 
                return self._problem.dJ(x)   
                
            def G(self, x): 
                return self._problem.G(x)   

            def dG(self, x): 
                return self._problem.dG(x)   

            def H(self, x):
                H = self._problem.H(x)
                if not isinstance(x, tuple):
                    x = (x,)
                for (xi,li,ui) in zip(x,l,u):
                    if not li is None:
                        H = np.hstack((H,li-xi))
                    if not ui is None:
                        H = np.hstack((H,xi-ui))
                return H
                
            def dH(self, x):
                dH = self._problem.dH(x)
                if isinstance(dH,list) and dH == []:
                    dH = None
                if not isinstance(x, tuple):
                    x= (x,)
                lengths = (0,)+tuple(len(xi) for xi in x)
                cumsums = np.cumsum(lengths)
                for (i,xi,li,ui) in zip(range(len(x)),x,l,u):
                    I = sp.eye(len(xi.flatten()), format="csc")
                    p = I.shape[0]
                    I = sp.hstack((sp.csc_matrix((p,cumsums[i])),
                                   I,
                                   sp.csc_matrix((p,cumsums[-1]-p-cumsums[i]))))
                    if not li is None:
                        dH = sp.vstack((dH,-I),format="csc")
                    if not ui is None:
                        dH = sp.vstack((dH, I), format="csc")
                if dH is None:
                    dH = []
                return dH

            def retract(self, x, dx):
                newx = self._problem.retract(x, dx)
                istuple = True
                if not isinstance(newx, tuple):
                    newx = (newx,)
                    istuple = False
                projected_newx =()
                for (xi, li, ui) in zip(newx, l, u):
                    if not li is None:
                        xi = np.maximum(xi, li)
                    if not ui is None:
                        xi = np.minimum(xi, ui)
                    projected_newx = projected_newx + (xi,)
                if not istuple:
                    return projected_newx[0]
                else:
                    return projected_newx
                
            def accept(self, params, results):
                results_copy = dict()   
                q = len(self._problem.H(results['x'][-1]))
                for key in results: 
                    if not key=='H':    
                        results_copy[key]=results[key]
                    else:   
                        results_copy[key] = [h[:q] for h in results['H']]   
                self._problem.accept(params, results_copy)
                
            def inner_product(self, x):
                return self._problem.inner_product(x)

        return BoundConstraintOptimizable
            
    return bound_constraints_optimizable_impl

