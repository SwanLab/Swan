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

import numpy as np
import scipy.sparse as sp
from .optimizable import Optimizable

class EuclideanOptimizable(Optimizable):
    r"""A subclass of :py:class:`~nullspace_optimizer.Optimizable`   
    for the optimization on Euclidean spaces, that is for solving   
    optimization problems of the form   

    .. math::   
        
       \begin{aligned} \min_{x\in \mathbb{R}^n} & \quad      J(x) \\
       s.t. & \left\{\begin{aligned}
        g_i(x) & =0  \text{ for all }0 \leqslant i \leqslant p-1,\\
        h_j(x) & \leqslant 0 \text{ for all }0\leqslant j \leqslant q-1,\end{aligned}\right.
        \end{aligned}

    where :math:`n` is the dimension of the design space,   
    :math:`p` and :math:`q` are the number of equality and inequality constraints.  
        
    On this case, the design variable is assumed to be a numpy array of length :math:`n`.
        
    Implementing a class that inherits this subclass should comply with the following structure:    

    .. code:: python
    
       from nullspace_optimizer import EuclideanOptimizable
    
       class MyEuclideanOptimizable(EuclideanOptimizable):
           # Initialization
           def x0(self):
               pass
    
           # Objective function
           def J(self, x):
               pass
    
           # Equality constraints
           def G(self, x):
               pass
    
           # Inequality constraints
           def H(self, x):
               pass
    
           # Derivative of the objective function
           def dJ(self, x):
               pass
    
           # Jacobian matrix of G
           def dG(self, x):
               pass
    
           # Jacobian matrix of H
           def dH(self, x):
               pass
    
           # Post processing every time a   
           # point on the optimization path is accepted
           def accept(self, params, results):
               pass

    It is not necessary to define the methods ``inner_product`` and ``retract`` 
    which are repsectively set automatically to the identity matrix 
    and to  
    ``retract(x,dx)=x + dx``.
    """

    def inner_product(self, x):
        dJ = np.asarray(self.dJ(x))
        if len(dJ.shape)>1:
            return sp.eye(dJ.shape[1], format = "csc")
        else:
            return sp.eye(len(dJ), format='csc')

    def retract(self, x, dx):
        return x+dx
