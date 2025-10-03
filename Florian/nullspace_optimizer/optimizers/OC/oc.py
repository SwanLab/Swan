import numpy as np
from nullspace_optimizer import Optimizable
from nullspace_optimizer.inout import display, tic, toc, display_iteration

from nullspace_optimizer.optimizers.utils import OptimizationResults, check_params

def oc_solve(problem : Optimizable, l, u, params = None, results = None):  
    r""" An Optimality Criteria solver for volume constrained compliance     
        minimization in topology optimization
        based on the script 
        `topopt_cholmod.py <https://www.topopt.mek.dtu.dk/-/media/subsites/topopt/apps/dokumenter-og-filer-til-apps/topopt_cholmod.py>`_.
            
        This solver applies to Optimization problems of the form    
                
        .. math::

           \begin{aligned}
               \min_{x\in \R^n}&  \quad J(x)\\
               \textrm{s.t.} & \left\{\begin{aligned}
            g(x)&=0,\\
            l  &\leqslant  x \leqslant u\\ 
                   \end{aligned}\right.
           \end{aligned}
            
        where :math:`x` is a numpy array and :math:`g` a single equality constraint. 
    """
    # Optimality Criteria Solver
    default_parameters = dict(maxit=2000,   
                              save_only_N_iterations = None,
                              move = 0.2,   
                              debug = 0)
    params = check_params(default_parameters, params)

    display("=================OC solver=================",
            color="magenta", attr="bold", level=1, debug=params['debug'])
    display("Params",
            color="magenta", attr="bold", level=5, debug=params['debug'])

    def oc(x,dc,dv,g):
        l1=0
        l2=1e9
        move=params['move']
        # reshape to perform vector operations
        xnew=np.zeros_like(x)

        while (l2-l1)/(l1+l2)>1e-3:
            lmid=0.5*(l2+l1)
            xnew[:]= np.maximum(l,np.maximum(x-move,np.minimum(u,np.minimum(x+move,x*np.sqrt(-dc/dv/lmid)))))
            gt=g+np.sum((dv*(xnew-x)))
            if gt>0 :
                l1=lmid
            else:
                l2=lmid
        return (xnew,gt)

    keys_group1 = ['it','J','G','H','x','g','change']
    keys_group2 = []
    results = OptimizationResults(keys_group1, keys_group2, results, params['save_only_N_iterations'])
        
    starting_values = results.initialize()
    if starting_values is None: 
        x = problem.x0()
        g = 0
        loop = 0    
        change = 1
    else:   
        x = starting_values['x']    
        loop = starting_values['it']
        change = starting_values['change']
        g = starting_values['g']

    xold=x.copy()

    J = problem.J(x)
    G = problem.G(x)[0]

    while change > 0.01 and loop < params['maxit']:    
        results.save('it',loop) 
        results.save('J',J)   
        results.save('G',[G])
        results.save('H',[])
        results.save('x',x.copy())
        results.save('g',g)
        results.save('change',change)
        display_iteration(loop, J, [G], [], x, level = 0, debug=params['debug'])
        problem.accept(params, results.implementation())

        loop = loop + 1
        dc = problem.dJ(x)  
        dv = problem.dG(x)
        xold[:]=x
        (x[:],g) = oc(x, dc, dv, g)

		# Compute the change by the inf. norm
        change=np.linalg.norm(x.flatten()-xold.flatten(),np.inf)

        J = problem.J(x)
        G = problem.G(x)[0]

    results.save('it',loop) 
    results.save('J',J)   
    results.save('G',[G])
    results.save('x',x)
    problem.accept(params, results.implementation())

    display('\n', -1, params['debug'])
    display('Optimization completed.', -1, params['debug'], color='blue')
    display_iteration(loop, J, [G], [], x, level = -1, debug=params['debug'])
    return results.implementation()        
