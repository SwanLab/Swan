from nullspace_optimizer import nlspace_solve, bound_constraints_optimizable, filtered_optimizable
from nullspace_optimizer import EuclideanOptimizable, memoize, minmax_optimizable
from nullspace_optimizer.examples.topopt_examples import ex00_compliance as ex00
from nullspace_optimizer.examples.basic_examples.utils import draw
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np

def init(nforces = 8): 
    global f
    n = np.ceil((ex00.nelx)/(nforces))  
    positions = np.array([k*n for k in range(nforces)]+[ex00.nelx],dtype=int)
    f = np.zeros((ex00.ndof,nforces))
    for i in range(nforces):    
        f[ex00.dofs[positions[i]:positions[i+1],0,1],i]= -1.0/ex00.nelx
    ex00.fixed = np.union1d(ex00.dofs[0,:,0],ex00.dofs[-1,-1,:])
    
@minmax_optimizable #Decorator for multiobjective optimization
@bound_constraints_optimizable(l=0,u=1)
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
        (obj,dc) = ex00.solve_state(x, f)
        return obj
        
    def dJ(self, x):
        (obj,dc)=ex00.solve_state(x, f)
        return dc

    def G(self,x):
        return [np.sum(x)/(ex00.nelx*ex00.nely)-ex00.volfrac]

    def dG(self,x):
        dv = np.ones(ex00.nelx*ex00.nely)/(ex00.nelx*ex00.nely)
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
    ex00.init(nelx=180,nely=90)
    init()
    case = TO_problem(plot=True)
    results = nlspace_solve(case,
                            {"dt": 0.1, 
                             'maxtrials':3,
                             'itnormalisation': 50,
                             'save_only_N_iterations':1,    
                             'save_only_Q_constraints':5,   
                             'qp_solver': 'qpalm',  
                             'tol_qp': 1e-8,
                             'maxit':150})
    plt.figure()
    draw.drawMuls(results)
    plt.figure()
    draw.drawJ(results)
    input("Press any key")


