# This file is part of nullspace_optimizer. 
#
# It is a modification of the 2018-2019 version initially developed with the    
# copyright of CNRS, Ecole Polytechnique and Safran.
#
# This version is maintained by the research team of Prof. Florian Feppon 2023
# https://people.cs.kuleuven.be/~florian.feppon/software.html
#
# nullspace_optimizer is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# nullspace_optimizer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# A copy of the GNU General Public License is included below.
# For further information, see <http://www.gnu.org/licenses/>.
import scipy.sparse as sp
from .optimizable import Optimizable
    

def filtered_optimizable(filter, diff_filter):  
    r"""A decorator for an Optimizable class for filtering optimization variables.   
        
    .. code:: python    
        
       def filter(x): 
            # Perform some operation on x   
            
       def diff_filter(x,v):    
           # returns v @ dfilter(x) where dfilter(x) is the Jacobian function   
           # of the filter

       @filtered_optimizable(filter, diff_filter)
       class problem(Optimizable):  
            def J(self, x): 
                #return something   
                    
            def G(self, x): 
                # return something  
                    
            #...

    Then, a problem of the form     

    .. math::
    
       \begin{aligned}
           \min_{x\in \R^n}&  \quad J(x)\\
           \textrm{s.t.} & \left\{\begin{aligned}
        g_i(x)&=0, \text{ for all } 1\leqslant i\leqslant p,\\
        h_j(x)  &\leqslant  0, \text{ for all }1\leqslant j \leqslant q,\\ 
               \end{aligned}\right.
       \end{aligned}

    is converted into 

    .. math::
    
       \begin{aligned}
           \min_{x\in \R^n}&  \quad J(\texttt{filter}(x))\\
           \textrm{s.t.} & \left\{\begin{aligned}
        g_i(\texttt{filter}(x))&=0, \text{ for all } 1\leqslant i\leqslant p,\\
        h_j(\texttt{filter}(x))  &\leqslant  0, \text{ for all }1\leqslant j \leqslant q,\\ 
               \end{aligned}\right.
       \end{aligned}

    Several decorators can be appended to implement successive composition of filters.   
        
    .. code:: python    

       @filtered_optimizable(filter_2, diff_filter_2)
       @filtered_optimizable(filter_1, diff_filter_1)
       class problem(Optimizable):  
       #...
        
    In that case, ``J(x)``, ``G(x)``
    will be replaced by ``J(filter_1(filter_2(x)))``, ``G(filter_1(filter_2(x)))``, etc...
        
    :param filter: the filtering function.
    :param diff_filter: the transpose derivative operator of the filtering function.  
                        ``diff_filter(x,v)`` sould be equal to ``v @ dfilter(x)``   
                        where ``dfilter`` is the Jacobian matrix of ``filter`` at ``x``, in 
                        order to match the chain rule   

                        .. code:: console

                           d[J(filter(x))] @ dx =dJ(filter(x)) @ dfilter(x) @ dx
                        
    """
    def filtered_optimizable_impl(optimizable): 
        class FilteredOptimizable(Optimizable):
            _undecorated = optimizable

            def __init__(self, *args, **kwargs): 
                self._problem = optimizable(*args, **kwargs) 
                self._xsave = []
                    
            def x0(self):   
                return self._problem.x0()

            def J(self, x): 
                return self._problem.J(filter(x))
                
            def G(self, x): 
                return self._problem.G(filter(x))   
                
            def H(self, x): 
                return self._problem.H(filter(x))

            def dJ(self, x):    
                dJ = self._problem.dJ(filter(x))  
                return diff_filter(x, dJ)

            def dG(self, x):    
                dG = self._problem.dG(filter(x)) 
                if dG is None or sp.csc_matrix(dG).shape[0]*sp.csc_matrix(dG).shape[1] == 0:
                    return []
                else:
                    return diff_filter(x,dG)

            def dH(self, x):    
                dH = self._problem.dH(filter(x)) 
                if dH is None or sp.csc_matrix(dH).shape[0]*sp.csc_matrix(dH).shape[1] == 0:
                    return []
                else:
                    return diff_filter(x,dH)
                
            def retract(self, x, dx):   
                return self._problem.retract(x, dx)
                
            def inner_product(self, x): 
                return self._problem.inner_product(filter(x))
                
            def accept(self, params, results):    
                results_copy = dict()
                self._xsave.append(filter(results['x'][-1]))
                if params['save_only_N_iterations'] and results['it'][-1]>=params['save_only_N_iterations']:    
                    self._xsave[-1-params['save_only_N_iterations']]=None

                for key in results: 
                    if not key=='x':    
                        results_copy[key] = results[key]
                    else:   
                        results_copy[key] = self._xsave
                self._problem.accept(params, results_copy)

        return FilteredOptimizable

    return filtered_optimizable_impl

