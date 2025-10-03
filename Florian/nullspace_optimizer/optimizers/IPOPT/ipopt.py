from ...optimizable import Optimizable  
from ...optimizers.utils import OptimizationResults, check_params
from ...inout import display_iteration, display
import numpy as np
from cyipopt import Problem
import scipy.sparse as sp
    
ipopt_lower_inf = -10**20
ipopt_upper_inf = 10**20

class IpoptProblem():    
    def __init__(self, problem: Optimizable, params, results : OptimizationResults): 
        self._problem = problem
        self.results = results  
        self.params = params
        
    def objective(self, x):     
        self.results.save('x',x)
        J = self._problem.J(x)  
        self.results.save('J',J)
        return self._problem.J(x)
        
    def gradient(self, x):  
        return self._problem.dJ(x)  
        
    def constraints(self, x):       
        G = self._problem.G(x)  
        H = self._problem.H(x)
        self.results.save('G',G)
        self.results.save('H',H)
        return np.hstack((G,H)) 
        
    def jacobian(self, x):  
        dG = self._problem.dG(x)
        dH = self._problem.dH(x)
        if hasattr(dG,'todense'): 
            dG = np.asarray(dG.todense())
        if hasattr(dH,'todense'): 
            dH = np.asarray(dH.todense())
        if dG is None or dG==[]:  
            return np.asarray(dH) 
        if dH is None or dH==[]:  
            return np.asarray(dG) 
        return np.vstack((dG,dH))   
        
    def intermediate(self, alg_mod, iter_count, obj_value, inf_pr, inf_du, mu,
                     d_norm, regularization_size, alpha_du, alpha_pr,
                     ls_trials):
        """Prints information at every Ipopt iteration."""
        self.results.save('it',iter_count)  
        for key in self.results.implementation().keys():    
            self.results.implementation()[key]=self.results[key][:iter_count+1]
        self._problem.accept(self.params,self.results.implementation())
        display_iteration(iter_count, self.results['J'][-1],    
                          self.results['G'][-1] if self.results['G'] else [],    
                          self.results['H'][-1] if self.results['H'] else [],    
                          self.results['x'][-1],  level = 0, debug= self.params['debug'])

def ipopt_solve(problem: Optimizable, l, u, params=None, results=None):   
    default_parameters = dict(debug = 0,
                              save_only_N_iterations=None,  
                              save_only_Q_constraints=None, 
                              ipopt_options=dict(), 
                              ipopt_print_level=0,
                              maxit=3000)
        
    params = check_params(default_parameters, params)
    params['ipopt_options']['max_iter'] = params['maxit']
    params['ipopt_options']['print_level'] = max(params['ipopt_print_level'],params['debug'])

    display("=================IPOPT (cyipopt interface)====================",
            color="magenta", attr="bold", level=1, debug=params['debug'])
    display("Params",
            color="magenta", attr="bold", level=5, debug=params['debug'])
    for key, value in params.items():
        display(f"{key} : {value}",
                color="magenta", level=5, debug=params['debug'])

    group1 = ['x', 'J', 'G', 'H', 'it']
    group2 = []
    abstract_results = OptimizationResults(group1, group2,
                                           results, params['save_only_N_iterations'],
                                           params['save_only_Q_constraints'])

    starting_values = abstract_results.initialize()
    if starting_values is None:
        x = problem.x0()
        it = 0
    else:
        x = starting_values['x']
        it = starting_values['it']
        display("Restart from iteration "+str(it),
                color="dark_olive_green_3b", attr="bold", level=0)
            
    if isinstance(l, list): 
        l = np.asarray(l)
    if isinstance(u, list):
        u = np.asarray(u)

    lb = l+np.zeros_like(x)
    ub = u+np.zeros_like(x) 

    G = problem.G(x)    
    H = problem.H(x)    
        
    p = len(G)  
    q = len(H)  
    cl = [0]*p+[ipopt_lower_inf]*q  
    cu = [0]*(p+q)
    nlp = Problem(n=len(x), 
                  m=len(cl),    
                  problem_obj=IpoptProblem(problem, params, abstract_results),
                  lb=lb,    
                  ub=ub,    
                  cl=cl,    
                  cu=cu)
    for key, value in params['ipopt_options'].items():  
        nlp.add_option(key, value)
    nlp.add_option('linear_solver', 'mumps')
        
    x, info = nlp.solve(x)
        
    return abstract_results.implementation()

