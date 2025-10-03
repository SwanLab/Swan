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

class EqualizedOptimizable(Optimizable):
    r"""An Optimizable object automatically converts all inequality constraints
    into equality constraints with the method of slack variables.

    **Usage**

    .. code:: python

       equalizedProblem = EqualizedOptimizable(problem)
        
    Supposing that ``problem`` corresponds to the optimization problem

    .. math::
    
       \begin{aligned}
           \min_{x\in \mathcal{X}}&  \quad J(x)\\
           \textrm{s.t.} & \left\{\begin{aligned}
        g_i(x)&=0, \text{ for all } 1\leqslant i\leqslant p,\\
        h_j(x)  &\leqslant  0 \text{ for all }1\leqslant j \leqslant q,\\ 
               \end{aligned}\right.
       \end{aligned}

    the ``equalizedProblem`` corresponds to the equalized program

    .. math::
    
       \begin{aligned}
           \min_{x\in \mathcal{X}}&  \quad J(x)\\
           \textrm{s.t.} & \left\{\begin{aligned}
        g_i(x)&=0, \text{ for all } 1\leqslant i\leqslant p,\\
        h_j(x)+\frac{1}{2}z_j^2  &= 0 \text{ for all }1\leqslant j \leqslant q,\\ 
               \end{aligned}\right.
       \end{aligned}

    Optimization points for the ``EqualizedOptimizable`` object are tuples of the
    form ``(x,z)`` with ``x`` an optimization point for `problem` and ``z`` a list of
    size :math:`q` the number of inequality constraints.

    The initial value for the slack variable is set to  

    .. math::   
        
       z_j(0) = \sqrt{2 |h_i(x(0))|} \text{ for } 1\leqslant j \leqslant q

    The inner product chosen on the z variable is the usual Euclidean inner
    product.  Derivatives with respect to the z variable are appended at the
    end of the derivative arrays with respect to x.
        
    .. note::   
        
       In practice, the equalization helps to find a central path going from    
       the initialization to the optimum in between the constraints. However,   
       equalization can be inefficient or require some tuning for problems with     
       a large number of constraints.
        
    :problem: An instance of :class:`Optimizable` implementing an arbitrary optimization    
              program.

    :param coeffs_z: some coefficients to equalize :math:`h_j(x)\leqslant 0` into   
                     :math:`h_j(x)+\frac{1}{2}\texttt{coeffs\_z}[j] z_j^2=0`.
    """

    def __init__(self, problem, coeffs_z = None):
        super().__init__()
        self.problem = problem
        self._results = dict()
        self._results['G'] = []
        self._results['H'] = []
        self.coeffs_z = coeffs_z    


    def x0(self):
        x0 = self.problem.x0()
        H = np.asarray(self.problem.H(x0))
        z = np.sqrt(2*np.abs(H))
        return (x0, z)

    def J(self, x):
        return self.problem.J(x[0])

    def G(self, x):
        oldG = self.problem.G(x[0])
        oldH = self.problem.H(x[0])
        Z = 0.5*x[1]**2
        if self.coeffs_z:   
            Z = Z*self.coeffs_z
        oldH = oldH + Z
        return np.hstack((oldG,oldH))

    def H(self, x):
        return []

    def dJ(self, x):
        return np.concatenate((self.problem.dJ(x[0]),
                               [0.0]*len(x[1])))

    def dG(self, x):
        old_dH = sp.csc_matrix(self.problem.dH(x[0]))
        z = x[1]
        if self.coeffs_z:   
            z = z*self.coeffs_z
        Z = sp.diags(z, format="csc")
        old_dH = sp.hstack((old_dH,Z))
        old_dG = self.problem.dG(x[0])  
        if old_dG is None or (isinstance(old_dG, list) and old_dG == []):
            return old_dH
        else:
            old_dG = sp.csc_matrix(old_dG)
        if old_dG.shape[0]*old_dG.shape[1] == 0:
            return old_dH   
        else:
            old_dG = sp.hstack((old_dG,sp.csc_matrix((old_dG.shape[0],len(x[1])))))
            return sp.vstack((old_dG,old_dH))

    def dH(self, x):
        return []

    def inner_product(self, x):
        Aold = self.problem.inner_product(x[0])
        A = sp.block_diag(
            (Aold, *(1,)*len(x[1])), format='csc')
        return A

    def retract(self, x, dx):
        retractedOld = self.problem.retract(
            x[0], dx[:(-len(x[1]))])
        retractedZi = x[1]+dx[-len(x[1]):]   
        return (retractedOld, retractedZi)

    def accept(self, params: dict, results: dict):
        for key in results: 
            if not key in {'x','G','H'}:
                self._results[key] = results[key]
            elif key == 'x':  
                self._results['x'] = [xi[0] if xi is not None    
                                     else None for xi in results['x']]
            elif key == 'G':
                # Save max 30 constraints
                self._results['G'].append(results['G'][-1][:min(30,-len(results['x'][-1][1]))])  
            elif key == 'H':  
                self._results['H'].append((results['G'][-1][-len(results['x'][-1][1]):]-0.5*results['x'][-1][1]**2)[:30])
        self.problem.accept(params, self._results)
