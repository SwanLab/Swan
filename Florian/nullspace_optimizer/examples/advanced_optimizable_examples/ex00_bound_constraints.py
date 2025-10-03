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

from nullspace_optimizer import EuclideanOptimizable, bound_constraints_optimizable
from nullspace_optimizer import nlspace_solve
import numpy as np
import matplotlib.pyplot as plt

#def checkx(x,n,title):  
#    print("Checking "+title)
#    N = int(np.sqrt(n))
#    x = np.reshape(x, (N, N))
#
#    cols = [x for x in reversed(range(N))]
#    print(x-x[cols,:])
#    assert np.allclose(x,x[cols,:])
#    assert np.allclose(x,x[:,cols])

@bound_constraints_optimizable(l=0,u=1)
class Problem(EuclideanOptimizable):
    def __init__(self, N, plot=None):
        self.N = N
        self.xi = np.linspace(0, 1, self.N)
        self.dx = self.xi[1]-self.xi[0]
        (self.X, self.Y) = np.meshgrid(self.xi, self.xi, indexing="ij")
        self.plot = plot
        if plot is None:
            self.plot = __name__=="__main__"

    def x0(self):
        x0 = np.abs(((self.X-0.5)*(self.Y-0.5)))    
        cols = [x for x in reversed(range(self.N))]
        halfN = int(self.N/2)
        x0[:,:halfN] = x0[:,cols[:halfN]]
        x0[:halfN,:] = x0[cols[:halfN],:]
        return x0.flatten()

    def J(self, x):
        #Maximize a discrete gradient energy norm to obtain a checkerboard
        x = np.reshape(x, (self.N, self.N))
        gradXsq = -0.5*(np.sum((x[1:, :]-x[:-1,:])**2) +
                        np.sum((x[:, 1:]-x[:, :-1])**2))
        return gradXsq

    def dJ(self, x):
        x = np.reshape(x, (self.N, self.N))
        gradJ = 0*x
        gradJ[1:, :] += x[1:, :]-x[:-1, :]
        gradJ[:-1, :] -= x[1:, :]-x[:-1, :]
        gradJ[:, 1:] += x[:, 1:]-x[:, :-1]
        gradJ[:, :-1] -= x[:, 1:]-x[:, :-1]
        return -gradJ.flatten()

    def accept(self, params, results):
        if self.plot:
            plt.figure(0)
            plt.figure(0).clear()
            plt.pcolor(self.X, self.Y, np.reshape(
                results['x'][-1], (self.N, self.N)), vmin=0, vmax=1, cmap="gray")
            #self.cbar = plt.colorbar()
            plt.figure(0).canvas.draw()
            plt.pause(0.1)
            

def run(N=10, plot=True, **options):
    #Bound constraints may be better dealt with 'normalize_tol':-1
    case = Problem(N, plot=plot)
    params = {"dt": 0.3, 'alphaJ': 1, 'alphaC': 1, 
              'itnormalisation': 10,   
              'save_only_N_iterations':1,    
              'save_only_Q_constraints':5,}  
    params.update(**options)
    results = nlspace_solve(case, params)
    return results  

def main():
    from nullspace_optimizer.examples.basic_examples.utils import draw
    results = run()
    plt.ion()

    plt.figure()
    draw.drawMuls(results)

    plt.figure()
    draw.drawJ(results)

    input("Press any key")

if __name__=="__main__":
    main()
