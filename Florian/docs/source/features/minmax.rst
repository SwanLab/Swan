Robust optimization 
-------------------
    
The ``nullspace_optimizer`` packages includes a decorator   
:func:`~nullspace_optimizer.minmax_optimizable` for conveniently implementing   
the robust min/max formulation of multiobjective optimization problems. 
    
It allows to implement the optimization problem     
    
.. math::

   \begin{aligned}
       \min_{x\in \mathcal{X}} \max_{1\leqslant i \leqslant r}& \quad J_i(x)\\
       \textrm{s.t.} & \left\{\begin{aligned}
    g_i(x)&=0, \text{ for all } 1\leqslant i\leqslant p,\\
    h_j(x)  &\leqslant  0, \text{ for all }1\leqslant j \leqslant q,\\ 
           \end{aligned}\right.
   \end{aligned}
    
by specifying a multivariate objective function in the Python implementation as
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
            
            

