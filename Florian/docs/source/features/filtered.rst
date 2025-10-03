Filtering of the design variable
--------------------------------
    
In some applications such as density based topology optimization, it can be useful to   
apply some transformation to the design variable. The   
decorator :func:`~nullspace_optimizer.filtered_optimizable` allows to automate this process, it     
converts an optimization problem 
    
.. math::

   \begin{aligned}
       \min_{x\in \R^n}&  \quad J(x)\\
       \textrm{s.t.} & \left\{\begin{aligned}
    g_i(x)&=0, \text{ for all } 1\leqslant i\leqslant p,\\
    h_j(x)  &\leqslant  0, \text{ for all }1\leqslant j \leqslant q,\\ 
           \end{aligned}\right.
   \end{aligned}
    
into the problem

.. math::

   \begin{aligned}
       \min_{x\in \R^n}&  \quad J(\rho(x))\\
       \textrm{s.t.} & \left\{\begin{aligned}
    g_i(\rho(x))&=0, \text{ for all } 1\leqslant i\leqslant p,\\
    h_j(\rho(x))  &\leqslant  0, \text{ for all }1\leqslant j \leqslant q,\\ 
           \end{aligned}\right.
   \end{aligned}
    
where :math:`\rho\,:\,\R^n\to\R^n` is a transformation of the design variable.  
    
The decorator is used as follows: 

.. code:: python    
    
   from nullspace_optimizer import EuclideanOptimizable, filtered_optimizable
        
   # Filter the design variable 
   # The optimizer will call J(filter(x)), G(filter(x)), etc.
   @filtered_optimizable(filter, diff_filter)
   class Problem(EuclideanOptimizable):  
        def x0(self):   
            #x0 
        #...
            
        def J(self, x): 
        # J 
            
        # etc...
    
    
The decorator takes two arguments:  
    
* ``filter``: the transformation map    
* ``diff_filter``: the derivative of the transformation map. It should be implemented such  
  that ``diff_filter(x,v)=v @ dfilter(x)`` where  ``dfilter(x)`` is the Jacobian matrix of the filter.
  This is required so that the decorator :func:`~nullspace_optimizer.filtered_optimizable`  
  can update consistently   
  the sensitivities ``DJ``, ``DG`` and ``DH`` according to the
  chain rule  
    
  .. math::

     \newcommand{\D}{\mathrm{D}}
     \frac{\D}{\D \texttt{x}}[\texttt{J}(\rho(\texttt{x}))]\cdot \texttt{dx} =
     \frac{\D \texttt{J}}{\D \texttt{x}}(\rho(\texttt{x})) \cdot \frac{\D \rho(\texttt{x})}{\D
     \texttt{x}}\cdot \texttt{dx}.
     

