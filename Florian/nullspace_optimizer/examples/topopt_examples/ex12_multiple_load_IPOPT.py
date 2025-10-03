from nullspace_optimizer import nlspace_solve, bound_constraints_optimizable, filtered_optimizable
from nullspace_optimizer import EuclideanOptimizable, memoize, minmax_optimizable, tuple_to_array_optimizable
from nullspace_optimizer.examples.topopt_examples import ex00_compliance as ex00
from nullspace_optimizer.optimizers.IPOPT import ipopt_solve, ipopt_lower_inf, ipopt_upper_inf
from nullspace_optimizer.examples.basic_examples.utils import draw
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import nullspace_optimizer.examples.topopt_examples.ex10_multiple_load as ex10
    
@tuple_to_array_optimizable
@minmax_optimizable
@filtered_optimizable(ex00.filter, ex00.diff_filter)
class TO_problem(EuclideanOptimizable):
    def __init__(self, plot=None):
        self.plot = plot
        if self.plot is None:    
            self.plot = __name__=="__main__"
        if self.plot:
            # Initialize plot and plot the initial design
            plt.ion() # Ensure that redrawing is possible
            self.fig,self.ax = plt.subplots()
            


    def x0(self):
        return ex00.volfrac * np.ones(ex00.nely*ex00.nelx,dtype=float)

    def J(self, x):
        (obj,dc) = ex00.solve_state(x, ex10.f)
        return obj
        
    def dJ(self, x):
        (obj,dc)=ex00.solve_state(x, ex10.f)
        return dc

    def H(self,x):
        return [np.sum(x)/(ex00.nelx*ex00.nely)-ex00.volfrac]

    def dH(self,x):
        dv = np.ones((1,ex00.nelx*ex00.nely))/(ex00.nelx*ex00.nely)
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
                

if __name__=="__main__":
    ex00.init(nelx=90,nely=30)
    ex10.init()
    case = TO_problem()
    l = [0]*ex00.nelx*ex00.nely + [ipopt_lower_inf]    
    u = [1]*ex00.nelx*ex00.nely + [ipopt_upper_inf]
    results = ipopt_solve(case, l, u,
                            {'save_only_Q_constraints':5,       
                             'maxit':200})

    plt.figure()
    draw.drawJ(results)
    input("Press any key")
