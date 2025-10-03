Memoization 
-----------
    
The ``nullspace_optimizer`` package includes the decorator  
:func:`~nullspace_optimizer.memoize`  for the memoization of functions. 
This allows to avoid multiple calls to an expensive function that needs to  
be used by several methods of an optimization problem, such as multiple constraints  
or the objective function and a constraint. 
    
**Usage**
    
Suppose we have ``func`` is an expensive function   
that returns a dictionary ``res`` with the value   
``res['J']`` of the objective function     
and a single constraint ``res['G']``. This function can be memoized as follows:
    
.. code:: python
    
   from nullspace_optimizer import memoize

   # enable memoization
   @memoize() 
   def func(x): 
       # An expensive function  
       # ....
       return res
        
Then a stored value will be returned if ``func`` is called several time     
with the same argument ``x``. This can be useful if ``func`` is needed in multiple  
methods of an :class:`Optimizable` object, such as ``J``, ``G`` or ``H``:  
    
.. code:: python    
    
   from nullspace_optimizer import Optimizable, nlspace_solve
    
   class Problem(Optimizable):  
    
       def x0(self):    
       # Initialization...  
        
       def J(self, x):  
          # Some expensive computation
          res = func(x) 
          return res['J'] 
            
       def G(self, x):  
          # Some expensive computation
          res = func(x)
          return res['G']
            
            
    # Solve the optimization problem    
    # Memoization will be used when calling J and G, func will be called    
    # only once at every iteration.
    results = nlspace_solve(Problem())


