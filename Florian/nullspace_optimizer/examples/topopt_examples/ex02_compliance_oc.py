# Copyright 2018-2019 CNRS, Ecole Polytechnique and Safran.
#
# This file is part of nullspace_optimizer.
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

# This test case is an adaptation of topopt_cholmod.py  from
# A 200 LINE TOPOLOGY OPTIMIZATION CODE BY NIELS AAGE
# AND VILLADS EGEDE JOHANSEN, JANUARY 2013
import nullspace_optimizer.examples.topopt_examples.ex00_compliance as ex00
from nullspace_optimizer.examples.topopt_examples.ex00_compliance import filtered_optimizable,\
                        filter, diff_filter, solve_state
from nullspace_optimizer.optimizers.OC import oc_solve
from nullspace_optimizer import EuclideanOptimizable
import matplotlib.pyplot as plt
from nullspace_optimizer.examples.basic_examples.utils import draw
import numpy as np
from matplotlib import colors

@filtered_optimizable(filter, diff_filter)
class TO_problem(EuclideanOptimizable):
    def __init__(self, plot=None):
        self.plot = plot
        if self.plot is None:    
            self.plot = __name__=="__main__"
        if self.plot:
            # Initialize plot and plot the initial design
            plt.ion() # Ensure that redrawing is possible
            self.fig,self.ax = plt.subplots()
        self.ex00 = ex00

    def x0(self):
        return ex00.volfrac * np.ones(ex00.nely*ex00.nelx,dtype=float)

    def J(self, x):
        (obj,dc)=solve_state(x)
        return obj[0]
        
    def dJ(self, x):
        (obj,dc)=solve_state(x)
        return dc[0,:]

    def G(self,x):
        return [np.sum(x)/(ex00.nelx*ex00.nely*ex00.volfrac)-1]

    def dG(self,x):
        dv = np.ones(ex00.nelx*ex00.nely)/(ex00.nelx*ex00.nely*ex00.volfrac)
        return dv

    def accept(self, params, results):
        if self.plot:
            x = results['x'][-1]

            if not hasattr(self,'im'):
                self.im = self.ax.imshow(-x.reshape((ex00.nelx,ex00.nely)).T, cmap='gray',\
                                    interpolation='none',norm=colors.Normalize(vmin=-1,vmax=0))
                self.ax.axis('off')
                self.fig.show()
            else:
                self.im.set_array(-x.reshape((ex00.nelx,ex00.nely)).T)
                self.fig.canvas.draw()
            plt.pause(0.01)

if __name__ == "__main__":
    ex00.init()
    case=TO_problem(plot=True)
    results = oc_solve(case, 0, 1, {'maxit':200,
                                    'save_only_N_iterations':1,    
                                    'move': 0.2})
    plt.figure()
    draw.drawMuls(results)
    plt.figure()
    draw.drawJ(results)
    input("Press any key")
