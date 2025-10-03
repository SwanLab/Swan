import nullspace_optimizer.examples.topopt_examples.ex07_heat_nlspace as ex07   
from nullspace_optimizer.examples.topopt_examples.ex07_heat_nlspace  \
    import filter, diff_filter, init, solve_state
from nullspace_optimizer import EuclideanOptimizable, filtered_optimizable  
from nullspace_optimizer.optimizers.OC import oc_solve
import matplotlib.pyplot as plt
import numpy as np
from pymedit import P0Function, P1Function
    
@filtered_optimizable(filter, diff_filter)
class Heat_TO(EuclideanOptimizable):    
    def __init__(self, plot=None): 
        self.plot = plot
        if self.plot is None:    
            self.plot = __name__=="__main__"
            
        if self.plot:   
            plt.ion()   
            self.fig, self.ax = plt.subplots(1,2)
            self.fig.set_size_inches(12,6)
            self.ax = self.ax.flatten()
            self.cbars = [None]*2
        
    def x0(self):   
        return ex07.volfrac * np.ones(ex07.Th.nt,dtype=float)
        
    def J(self, rho): 
        return solve_state(rho)['J']
        
    def G(self, rho): 
        return [solve_state(rho)['H']]
        
    def dJ(self, rho):    
        return solve_state(rho)['dJ']
        
    def dG(self, rho):    
        return solve_state(rho)['dH']
        
    def accept(self, params, results):    
        if self.plot:
            self.ax[0].clear()  
                
            rho = results['x'][-1]
            for i in range(2):  
                self.ax[i].clear()  
                if self.cbars[i]:
                    self.cbars[i].remove()
            _, _, self.cbars[0] = P0Function(ex07.Th, rho).plot(fig=self.fig, ax=self.ax[0],cmap="gray_r", 
                                                         vmin=0, vmax=1)
            _, _, self.cbars[1] = P1Function(ex07.Th, solve_state(rho)['T[]']).plot(fig=self.fig, cmap="turbo",    
                                                                             ax=self.ax[1], 
                                                       type_plot='tricontourf')
            self.ax[0].title.set_text('rho')
            self.ax[1].title.set_text('T')
            self.fig.suptitle('Iteration '+str(results['it'][-1]))
            # savefig(fig_dir+"/"+it+".png", fig=self.fig, level=1, sizefactor=1)
            plt.pause(0.1)

if __name__=="__main__":    
    init(100)
    plt.ion()

    case = Heat_TO(plot=True)
    params = dict(move=0.3,   
                  maxit=150,    
                  save_only_N_iterations=3)
    oc_solve(case, 0, 1, params)
