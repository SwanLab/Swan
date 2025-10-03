Bound constraints
-----------------
    
The ``nullspace_optimizer`` package includes the decorator  
:func:`~nullspace_optimizer.bound_constraints_optimizable`  
for automating the implementation of bound constraints. It converts  
an optimization problem of the form     
    
.. math::

   \begin{aligned}
       \min_{x\in \R^n}&  \quad J(x)\\
       \textrm{s.t.} & \left\{\begin{aligned}
    g_i(x)&=0, \text{ for all } 1\leqslant i\leqslant p,\\
    h_j(x)  &\leqslant  0, \text{ for all }1\leqslant j \leqslant q,\\ 
           \end{aligned}\right.
   \end{aligned}
    
into 

.. math::

   \begin{aligned}
       \min_{x\in \R^n}&  \quad J(x)\\
       \textrm{s.t.} & \left\{\begin{aligned}
    g_i(x)&=0, \text{ for all } 1\leqslant i\leqslant p,\\
    h_j(x)  &\leqslant  0, \text{ for all }1\leqslant j \leqslant q,\\ 
    l& \leqslant x \leqslant u,
           \end{aligned}\right.
   \end{aligned}
    
where :math:`l` and :math:`u` are lower and upper bound vectors of dimension :math:`n`.
The decorator is used as follows:   
    
.. code:: python    
    
   from nullspace_optimizer import EuclideanOptimizable, bound_constraints_optimizable
        
   # Append the bound constraints l <= x <= u on the    
   # design variable x
   @bound_constraints_optimizable(l,u)
   class Problem(EuclideanOptimizable):  
        def x0(self):   
            #x0 
        #...
            
        def J(self, x): 
        # J 
            
        # etc...
            
The decorator :func:`~nullspace_optimizer.bound_constraints_optimizable`   
appends the bound constraints :math:`\texttt{l}\leqslant \texttt{x}\leqslant \texttt{u}` 
to the optimization problem ``Problem``, supposing that the design variable     
``x``, the lower bound ``l`` and the upper bound ``u`` are numpy arrays of the same size.
    
The design variable can also be a tuple of numpy array. See the full    
documentation of the decorator :func:`~nullspace_optimizer.bound_constraints_optimizable`.
