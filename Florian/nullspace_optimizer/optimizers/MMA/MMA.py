import numpy as np
from ...optimizable import Optimizable
from ...inout import display, tic, toc, display_iteration

from ..utils import OptimizationResults, check_params
from .mma_lib import mmasub, subsolv
from . import mma_lib 

def reshape(x): 
    try:
        return x.flatten()[np.newaxis].T
    except:
        import ipdb 
        ipdb.set_trace()

def checkx(x,n,title):  
    print("Checking "+title)
    N = int(np.sqrt(n))
    x = np.reshape(x, (N, N))

    cols = [x for x in reversed(range(N))]
    print(x-x[cols,:])
    assert np.allclose(x,x[cols,:])
    assert np.allclose(x,x[:,cols])

def mma_solve(problem: Optimizable, l, u, params=None, results=None):
    default_parameters = dict(maxit=2000,
                              tol=1e-3,
                              save_only_N_iterations=None,
                              save_only_Q_constraints=None, 
                              a0 = 1, 
                              a = None, 
                              c = None, 
                              d = None, 
                              fix = False,
                              move = 0.1,   
                              tight_move = False,
                              debug = 0)
    params = check_params(default_parameters, params)

    display("=================MMA=================",
            color="magenta", attr="bold", level=1, debug=params['debug'])
    display("Params",
            color="magenta", attr="bold", level=5, debug=params['debug'])
    for key, value in params.items():
        display(f"{key} : {value}",
                color="magenta", level=5, debug=params['debug'])

    mma_lib.FIX = params['fix']


    H = problem.H(problem.x0())
    m = len(H)
    for key in ['a','d']:
        if params[key] is None: 
            params[key] = np.zeros((m,1))
    if params['c'] is None: 
        params['c'] = np.ones((m,1))*10000
    else:   
        params['c'] = params['c']+np.zeros((m,1))
                
    if params['save_only_N_iterations'] == 1:    
        raise Exception("Error: params['save_only_N_iterations'] cannot be equal to 1"   
                        " for MMA to run correctly.")

    group1 = ['x', 'J', 'H', 'it', 'change', 'low', 'upp','muls']
    group2 = []
    abstract_results = OptimizationResults(group1, group2,
                                           results, params['save_only_N_iterations'],
                                           params['save_only_Q_constraints'])

    starting_values = abstract_results.initialize()
    if starting_values is None:
        x = problem.x0()
        n = len(x)
        it = 0
        change = 1
        low = np.ones((n,1))
        upp = np.ones((n,1))
        lam = None  
    else:
        x = starting_values['x']
        it = starting_values['it']
        change = starting_values['change']
        low = starting_values['low']
        upp = starting_values['upp']
        lam = starting_values['muls']
        display("Restart from iteration "+str(it),
                color="dark_olive_green_3b", attr="bold", level=0)

    J = problem.J(x)    
    H = problem.H(x)
        
    Jscale = J/10

    xval = reshape(x)
    if isinstance(l, list): 
        l = np.asarray(l)[:,np.newaxis]
    if isinstance(u, list):
        u = np.asarray(u)[:,np.newaxis]

    Xmin = l+np.zeros_like(xval)
    Xmax = u+np.zeros_like(xval)

    while change > params['tol'] and it < params['maxit']:
        f0val = J / Jscale
        df0dx = problem.dJ(x)[np.newaxis].T / Jscale
        fval = np.array(H)
        dfdx = problem.dH(x)
        if hasattr(dfdx,'todense'): 
            dfdx = np.asarray(dfdx.todense())
        if m==1:    
            dfdx = dfdx[np.newaxis]
        abstract_results.save('J', J)
        abstract_results.save('H', H)
        abstract_results.save('it', it)
        abstract_results.save('change', change)
        abstract_results.save('x', x)
        abstract_results.save('muls', lam)
        abstract_results.save('low', low)
        abstract_results.save('upp', upp)
        display_iteration(it, J, [], fval, x, level = 0,  debug = params['debug'])
        problem.accept(params, abstract_results.implementation())

        xval = reshape(x)
        if params['tight_move']:
            xmin = np.maximum(Xmin,xval-params['move'])
            xmax = np.minimum(Xmax,xval+params['move'])
        else:
            xmin = Xmin 
            xmax = Xmax
        m = len(fval)
        n = len(xval)
        if it >= 1:
            xold1 = reshape(abstract_results['x'][-2])
        else:   
            xold1 = reshape(x)
        if it >= 2: 
            xold2 = abstract_results['x'][-3]
            xold2 = reshape(xold2)
        else:   
            xold2 = xold1
            
        xmma, ymma, zmma, lam, xsi, eta, mu, zet, s, low, upp = \
            mmasub(m, n, it+1, xval, xmin, xmax, xold1, xold2, f0val,
                   df0dx, fval[:,np.newaxis], dfdx, low, upp,     
                   params['a0'], params['a'], params['c'], params['d'], 
                   params['move'])
            
        x = xmma.flatten()
        change = np.linalg.norm(xmma.flatten()-xold1.flatten(),np.inf)
        it = it+1
        J = problem.J(x)    
        H = problem.H(x)

    abstract_results.save('J', J)
    abstract_results.save('H', H)
    abstract_results.save('it', it)
    abstract_results.save('change', change)
    abstract_results.save('x', x)
    abstract_results.save('muls', lam)
    problem.accept(params, abstract_results.implementation())

    display('\n', -1, params['debug'])
    display('Optimization completed.', -1, params['debug'], color='blue')
    display_iteration(it, J, [], H, x,  level = -1, debug= params['debug'], color='blue')
    return abstract_results.implementation()
            

