.. _symbolic_differentiation:    

Symbolic differentiation    
-------------------------
    
The ``nullspace_optimizer`` package includes the decorator  
:func:`~nullspace_optimizer.symbolic_optimizable`  
that enables to  implement optimization programs with     
symbolic objective and constraint functions.    
    
For instance, the optimization problem  

.. math::       
    
   \newcommand{\<}{\leq}
   \begin{aligned} \min_{(x_0,x_1)\in\mathbb{R}^2} & \quad x_1+0.3 x_0  \\
   s.t. &\quad  \left\{ \begin{aligned}  -x_1+\frac{1}{x_0} & \< 0\\
                                -(3-x_0-x_1) & \< 0 
                                \end{aligned}\right.
   \end{aligned}

can be implemented symbolically as follows: 
    
.. code:: python

   from nullspace_optimizer import EuclideanOptimizable, symbolic_optimizable

   @symbolic_optimizable
   class basicProblem(EuclideanOptimizable):
       def x0(self):
           return [1.5, 2.25]
   
       def J(self, x):
           return x[1]+0.3*x[0]
   
       def H(self, x):
           return [-x[1]+1.0/x[0], -(3-x[0]-x[1])]


Only the objective and constraint functions need to be specified,   
their derivatives being obtained with   
`Sympy <https://docs.sympy.org/latest/index.html>`_ symbolic differentiation.  
The decorator :func:`~nullspace_optimizer.symbolic_optimizable`  
is compatible with the use of arbitrary `Sympy <https://docs.sympy.org/latest/index.html>`_ 
symbolic operators.
