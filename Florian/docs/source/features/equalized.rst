Equalization of inequality constraints  
--------------------------------------
    

The ``nullspace_optimizer`` package includes the class  
:class:`~nullspace_optimizer.EqualizedOptimizable`  
that converts inequality constraints into equality constraints by adding slack variables.   
More precisely, it converts  
an optimization problem of the form     
    
.. math::

   \begin{aligned}
       \min_{x\in \mathcal{X}}&  \quad J(x)\\
       \textrm{s.t.} & \left\{\begin{aligned}
    g_i(x)&=0, \text{ for all } 1\leqslant i\leqslant p,\\
    h_j(x)  &\leqslant  0, \text{ for all }1\leqslant j \leqslant q,\\ 
           \end{aligned}\right.
   \end{aligned}
    
into 

.. math::

   \begin{aligned}
       \min_{x\in \mathcal{X}}&  \quad J(x)\\
       \textrm{s.t.} & \left\{\begin{aligned}
    g_i(x)&=0, \text{ for all } 1\leqslant i\leqslant p,\\
    h_j(x)+\frac{1}{2}z_j^2  &= 0, \text{ for all }1\leqslant j \leqslant q,\\ 
           \end{aligned}\right.
   \end{aligned}
    
The class is used as follows:
    
.. code:: python    
    
   from nullspace_optimizer import Optimizable, EqualizedOptimizable, \
                                   nlspace_solve
        
   class Problem(Optimizable):  
        def x0(self):   
            #x0 
        #...
            
        def J(self, x): 
        # J 
            
        # etc...
            
    # Equalized version of the optimization problem Problem
    equalized_problem = EqualizedOptimizable(Problem())
    # Solve the equalized version of the problem
    results = nlspace_solve(equalized_problem)
