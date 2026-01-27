## PACKAGES
import sys
sys.path.append("/home/joseantonio/Documentos/GitHub/Swan/Florian"), \

from nullspace_optimizer import EuclideanOptimizable,\
bound_constraints_optimizable, memoize, filtered_optimizable
import numpy as np
from nullspace_optimizer.optimizable import Optimizable
from nullspace_optimizer.inout import tic, toc   
from nullspace_optimizer import inout as io   
from nullspace_optimizer.optimizers.nullspace.utils import compute_norm, get_xiJ_xiC, get_tilde, pack_constraints
from nullspace_optimizer.optimizers.utils import OptimizationResults, check_params
from nullspace_optimizer.optimizers.nullspace import utils
import cvxopt
import cvxopt.cholmod
import scipy.sparse as sp
import matplotlib.pyplot as plt
from matplotlib import colors
from pyfreefem import FreeFemRunner
from pymedit import P1Function



def FunctionCase04(case,No,maxIter):

    ## FREEFEM PROBLEM DEFINITION
    path = "NonLinearNull/04_NullSpace_HJ_Compl_VolPer_Newton_OldDual/"

    exports = FreeFemRunner(path+"04_Mesh.edp").execute()
    Th = exports['Th']
    alpha = exports['alpha']
    beta = exports['beta']
    labelDir = exports['labelDir']
    labelNeu = exports['labelNeu']
    meshsiz = exports['meshsiz']
    lsLabel = 10
    rInner = 3

    @bound_constraints_optimizable()
    class TO_problem(EuclideanOptimizable):
        def __init__(self):
            self.volFrac = []
            self.ux = []
            self.uy = []
            self.Th2 = []
            self.nx = []
            self.ny = []
        def x0(self):
            runner = FreeFemRunner(path+"04_InitialGuess.edp")
            runner.import_variables(Th=Th)
            return runner.execute()['phi[]']

        def J(self, x):
            runner = FreeFemRunner(path+"04_Cost.edp")
            runner.import_variables(Th=Th,phiVal=x,labelDir=labelDir,labelNeu=labelNeu,
                                    AchiVal=self.volFrac)
            exports = runner.execute()
            self.ux = exports['ux[]']
            self.uy = exports['uy[]']
            return exports['J']

        def dJ(self, x):
            runner = FreeFemRunner(path+"04_CostGradient.edp")
            runner.import_variables(Th=Th,Th2=self.Th2,beta=beta,lsLab=lsLabel,
                                    uxVal=self.ux,uyVal=self.uy)
            return runner.execute()['g[]']

        def G(self, x):
            runner = FreeFemRunner(path+"04_ConstraintEq.edp")
            runner.import_variables(Th=Th,AchiVal=self.volFrac)
            return [runner.execute()['C']]

        def dG(self, x):
            runner = FreeFemRunner(path+"04_ConstraintEqGradient.edp")
            runner.import_variables(Th=Th,Th2=self.Th2,phiVal=x,beta=beta,
                                    lsLab=lsLabel,rInner=rInner)
            return runner.execute()['g[]']
        
        def H(self, x):
            runner = FreeFemRunner(path+"04_ConstraintIneq.edp")
            runner.import_variables(Th=Th,phiVal=x)
            return [runner.execute()['H']]

        def dH(self, x):
            runner = FreeFemRunner(path+"04_ConstraintIneqGradient.edp")
            runner.import_variables(Th=Th,Th2=self.Th2,phiVal=x,nxVal=self.nx,
                                    nyVal=self.ny,alpha=alpha,beta=beta,
                                    lsLab=lsLabel,rInner=rInner)
            return runner.execute()['g[]']

        def accept(self, params, results):
            # Plot the design at every iteration
            x = results["x"][-1]
            u = P1Function(Th,x<=0)
            fig = params['fig']
            ax = params['ax']
            ax.clear()
            u.plot(
            fig=fig,
            ax=ax,
            title="Design variable"
            )
            plt.pause(0.05)

    ## OPTIMIZATION PARAMETERS
    dTime = 0.001
    hmin = meshsiz
    elRadius = 10
    params = {"dt": dTime*hmin*elRadius,
            "itnormalisation": No,
            "save_only_N_iterations": 1,
            "save_only_Q_constraints": 5,
            "alphaJ": 1,
            "alphaC": 1,
            "maxit": maxIter,
            "CFL": 0.9}
    problem:Optimizable = TO_problem()




    ## SETTINGS MANAGEMENT
    default_parameters = dict(alphaJ=1,  
                                alphaC=1,     
                                maxit=4000,    
                                debug=0,   
                                normalisation_norm=np.inf, 
                                dual_norm = 1,
                                itnormalisation=1, 
                                tol_finite_diff=0.15,  
                                dt=0.1,    
                                K=0.1, 
                                tol_cg=1e-15,  
                                alphas=lambda x : None,    
                                maxtrials=3,
                                normalize_tol = -1,
                                qp_solver_options = dict(),    
                                qp_solver = 'osqp',    
                                show_progress_qp = False,  
                                method_xiC = "linear_system",
                                tol_qp = 1e-10,    
                                CFL = 0.9,
                                provide_gradients = False,
                                save_only_N_iterations = None,     
                                save_only_Q_constraints = None,   
                                start = None)

    results = None

    if not params is None and 'dt' in params: 
        default_parameters['tol'] = 1e-5*params.get('dt') 
    else:   
        default_parameters['tol'] = 1e-5*default_parameters['dt']
        
    params = check_params(default_parameters, params)

    io._display.debug = params['debug']
    if params['debug'] >= 10:   
        params['show_progress_qp'] = True

    if params['qp_solver'] == 'cvxopt':   
        params['qp_solver_options'].update({'show_progress':params['show_progress_qp'], 
                                'abstol': params['tol_qp']})
    elif params['qp_solver'] == 'osqp':   
        params['qp_solver_options'].update({'verbose':params['show_progress_qp'],   
                                'eps_abs': params['tol_qp']})
    elif params['qp_solver'] == 'qpalm':   
        params['qp_solver_options'].update({'verbose':params['show_progress_qp'],   
                                'eps_abs': params['tol_qp']})
    else:   
        raise Exception("Wrong qp solver provided: "+params['qp_solver']+". Available: cvxopt|osqp|qpalm")
    utils.tol_qp = params['tol_qp']
    utils.tol_cg = params['tol_cg']
        
    if not callable(params['alphas']):  
        alphas = params['alphas']
        params['alphas'] = lambda x : alphas
                                
    io.display("="*20+"Null Space Optimizer"+"="*20,
            color="magenta", attr="bold", level=1)
    io.display("Params",
            color="magenta", attr="bold", level=5)
    for key, value in params.items():
        io.display(f"{key} : {value}",
                color="magenta", level=5)

    if params['provide_gradients']: 
        def get_gradient_transpose(A, dJ, dC, tildeEps):    
            dJT = problem.dJT(x)    
            if p==0:   
                dGT = np.empty((len(dJ),0))
            else:
                dGT = problem.dGT(x)
            if q==0:
                dHT = np.empty((len(dJ),0))
            else:
                dHT = problem.dHT(x)
            dCT = np.hstack((dGT,dHT))
            return dJT, dCT
        utils.get_gradient_transpose = get_gradient_transpose 

    plt.ion()
    fig, ax = plt.subplots()
    params['fig'] = fig
    params['ax'] = ax

    # Load previous results
    group1 = ['x','J', 'G', 'H', 's', 'it','normxiJ_save','muls']
    group2 = ['normxiJ', 'eps', 'tolerance']
    abstract_results = OptimizationResults(group1, group2,  
                                        results, params['save_only_N_iterations'], 
                                        params['save_only_Q_constraints'],   
                                        params['start'])

    starting_values = abstract_results.initialize()
    if starting_values is None: 
        x = problem.x0()
        it = 0  
        s = 0   
        normxiJ_save = None     
        muls = None
    else:
        x = starting_values['x']    
        it = starting_values['it']
        s = starting_values['s']    
        normxiJ_save = starting_values['normxiJ_save']  
        muls = starting_values['muls']
        io.display("Restart from iteration "+str(it), color="dark_olive_green_3b", attr="bold", level=0)






    ## NULLSPACE ALGORITHM

    runner = FreeFemRunner(path+"04_VolFracComputer.edp")
    runner.import_variables(Th=Th,phiVal=x)
    problem._problem.volFrac = runner.execute()['Achi[]']

    J = problem.J(x)
    G = problem.G(x)
    H = problem.H(x)

    normdx = 1  # current value for x_{n+1}-x_n
    if muls is None:
        muls = np.zeros(len(G) + len(H))

    #Optimization loop
    while normdx > params['tol'] and it <= params['maxit']:
        abstract_results.save('it', it)
        abstract_results.save('J', J)
        abstract_results.save('G', G)
        abstract_results.save('H', H)
        abstract_results.save('x', x)
        abstract_results.save('s', s)
        abstract_results.save('muls', muls)
        abstract_results.save('normxiJ_save', normxiJ_save)
        tic()
        problem.accept(params, abstract_results.implementation())
        io.display("Accept took "+toc(), level=4)
        x = abstract_results['x'][-1]

        io.display_iteration(it, J, G, H, x, level = 0,  debug = params['debug'])

        runner = FreeFemRunner(path+"04_BoundaryRefinement.edp")
        runner.import_variables(Th=Th,phiVal=x,alpha=alpha,lsLab=lsLabel,rInner=rInner)
        exports = runner.execute()

        Th2Old = exports['Th2']
        nxOld = exports['nx[]']
        nyOld = exports['ny[]']
        
        problem._problem.Th2 = Th2Old
        problem._problem.nx = nxOld
        problem._problem.ny = nyOld

        A = problem.inner_product(x)

        dJ = problem.dJ(x)

        itj = 1
        xj = x

        dG = problem.dG(xj)
        dH = problem.dH(xj)

        J, G, H, dJ, dG, dH, C, dC, n, p, q = pack_constraints(J, G, H, dJ, dG, dH)

        prevmuls = muls.copy()

        # Null space direction xiJ and range space direction xiC
        tic()
        xiJ, xiC, eps, tildeEps, muls = get_xiJ_xiC(J, G, H, dJ, dG, dH, A,  
                                            h=params['K']*params['dt'],        
                                            alphas = params['alphas'](xj),
                                            dual_norm=params['dual_norm'],   
                                            qp_solver = params['qp_solver'],     
                                            qp_solver_options = params['qp_solver_options'], 
                                            method_xiC = params['method_xiC'])
            
        abstract_results.save('eps', eps)
        io.display(f"Lagrange multipliers: {muls[:10]}", 5, params['debug'], "magenta")

        normxiJ = compute_norm(xiJ, params['normalisation_norm'])
        abstract_results.save('normxiJ', normxiJ)
                
        # Rescalings 
        if it < params['itnormalisation'] or \
        (params['normalize_tol'] >= 0 and
            not np.all((muls[p:] > params['normalize_tol'])
                    == (prevmuls[p:] > params['normalize_tol']))):
            io.display(f"Normalisation of the xiJ direction, params['itnormalisation']={max(params['itnormalisation'], it + 1)}",
                     level=5)
            normxiJ_save = normxiJ
            AJ = (params['alphaJ'] * params['dt']) / (1e-9 + normxiJ)
        else:
            AJ = params['alphaJ']*params['dt'] / max(1e-9 + normxiJ, normxiJ_save)
        AC = min(params['CFL'], params['alphaC']*params['dt'] / max(compute_norm(xiC, params['normalisation_norm']), 1e-9))

        # Update step
        g = AJ*xiJ+AC*xiC
        runner = FreeFemRunner(path+"04_HJ.edp")
        runner.import_variables(Th=Th,Th2=problem._problem.Th2,gVal = g,phiVal=xj,nxVal=problem._problem.nx,
                                    nyVal=problem._problem.ny,beta=beta,lsLab=lsLabel,rInner=rInner,dTime=dTime)
        x1 = runner.execute()['phi[]']

        runner = FreeFemRunner(path+"04_BoundaryRefinement.edp")
        runner.import_variables(Th=Th,phiVal=x1,alpha=alpha,lsLab=lsLabel,rInner=rInner)
        exports = runner.execute()

        problem._problem.Th2 = exports['Th2']
        problem._problem.nx = exports['nx[]']
        problem._problem.ny = exports['ny[]']

        #if np.linalg.norm(G,ord=np.inf)<0.01 and max(H)<0.01:
            #maxItj = 1
            #params['alphaJ'] = 1
        dx = (x1-x)
        #if p>0:
        #    assert np.isclose(dC[:p,:] @ xiJ,0,atol=1e-15) 
        #if dC[tildeEps,:][p:,:].size > 0:
        #    print(np.min(dC[tildeEps,:][p:,:]@xiJ))

        # Tolerance bounds at which one can expect to meet the constraint
        abstract_results.save('tolerance',  
                            np.array(np.sum(abs(dC), 1))*compute_norm(dx,np.inf))

        # Loop with trials
        normdx = compute_norm(dx, 2)
        success = 0
        tilde = get_tilde(C, p)
        for k in range(params['maxtrials']):
            newx = problem.retract(x, (0.5**k)*dx)

            runner = FreeFemRunner(path+"04_VolFracComputer.edp")
            runner.import_variables(Th=Th,phiVal=newx)
            problem._problem.volFrac = runner.execute()['Achi[]']

            (newJ, newG, newH) = (problem.J(newx), problem.G(newx), problem.H(newx))
            newC = np.concatenate((newG, newH))
                
            # Finite difference check
            finite_diffJ = np.abs(
                newJ - J - dJ.dot((0.5**k)*dx))/(0.5**k*normdx+1e-10)
            finite_diffC = max(
                np.abs(newC - C - dC.dot((0.5**k)*dx))/(0.5**k*normdx+1e-10), default=0)
            if max(finite_diffJ, finite_diffC) > params['tol_finite_diff']:
                io.display("Warning, inaccurate finite differences, time step might be too large. "
                        f"finite_diffJ={finite_diffJ}, finite_diffC={finite_diffC}", 1, params['debug'], color="dark_orange_3a")
            if newJ > J and np.linalg.norm(newC[tilde],2) >= np.linalg.norm(C[tilde],2):
                io.display(f"Warning, newJ={newJ} > J={J} and normNewC={np.linalg.norm(newC[tilde],2)} > normC= {np.linalg.norm(C[tilde],2)} "
                        + f"-> Trial {k+1}", 0, params['debug'], color="red")
            else:
                success = 1
                break

        if not success:
            io.display(
                "All trials have failed, passing to the next iteration.", 0, params['debug'],
                color="red")

        x = newx
        it += 1
        (J, G, H) = (newJ, newG, newH)
        # Optimization path length
        s += (0.5**k)*np.linalg.norm(dx, 2)

    abstract_results.save('J', J)
    abstract_results.save('G', G)
    abstract_results.save('H', H)
    abstract_results.save('x', x)
    abstract_results.save('it', it)
    abstract_results.save('s', s)
    abstract_results.save('muls', muls)
    abstract_results.save('normxiJ_save', normxiJ_save)
    problem.accept(params, abstract_results.implementation())

    io.display('\n', -1, params['debug'])
    io.display('Optimization completed.', -1, params['debug'], color='blue')
    io.display_iteration(it, J, G, H, x,  level = -1, debug= params['debug'], color='blue')




    ## POSTPROCESS
    results = abstract_results.implementation()
    iter = results['it']
    Comp  = results['J']
    Vol  = results['G']
    Per = results['H']

    np.savez(path+"04_ResultCase"+str(case),
            xF=x,it=iter,c=Comp,v=Vol,p=Per)