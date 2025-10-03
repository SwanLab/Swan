import numpy as np
from ...optimizable import Optimizable
from ...inout import tic, toc   
from ... import inout as io   
from .utils import compute_norm, get_xiJ_xiC, get_tilde, pack_constraints
from ..utils import OptimizationResults, check_params
from . import utils


def nlspace_solve(problem: Optimizable, params=None, results=None):
    r"""
    Null Space Optimizer solver for the nonlinear constrained optimization problem

    .. math::   
       :label: opt_pb
        
       \begin{aligned} \min_{x\in \mathcal{X}} & \quad      J(x) \\
       s.t. & \left\{\begin{aligned}
        g_i(x) & =0  \text{ for all }0 \leqslant i \leqslant p-1,\\
        h_j(x) & \leqslant 0 \text{ for all }0\leqslant j \leqslant q-1,\end{aligned}\right.
        \end{aligned}

    where :math:`p` and :math:`q` are the number of equality and inequality constraints. The    
    Null Space Optimizer solves the optimization problem :eq:`opt_pb` by solving the    
    Ordinary Differential Equation (so-called ‘’null-space gradient flow’’)

    .. math::
       :label: nlspace_ode
    
        \dot{x}(t) =-\alpha_J(t) \xi_J(x(t))-\alpha_C(t)\xi_C(x(t))
        
    with adaptive coefficients :math:`\alpha_J(t)` and :math:`\alpha_C(t)`.

    :param problem: An :py:class:`~nullspace_optimizer.Optimizable` 
                    object implementing the optimization problem. 

    :param params: (optional) a dictionary containing algorithm parameters
                    (see below).

    :param results: (optional) a previous output of the :py:class:`~nullspace_optimizer.nlspace_solve`  
                    function. The optimization will restart from the last input of
                    the dictionary ``results['x'][-1]`` and keep going as if there were   
                    no interruption.

    :return results: A dictionary of results containing the following entries:

                     - ``results['it']``: iteration numbers ``[0, 1, ..., n]``. 
                     - ``results['x']``: optimization path ``[x0,x1,...,xn]``
                     - ``results['J']``: values ``[J(x_0),...,J(x_n)]`` of the objective function along the path
                     - ``results['G']``: equality constraint values
                       ``[G(x0),...,G(xn)]``
                     - ``results['H']``: inequality constraints values
                       ``[H(x_0),...,H(x_n)]``
                     - ``results['muls']``: lagrange multiplier values at every iterations
                       ``[mu(x0),...,mu(xn)]``. Each ``mu(xi)`` is a vector or list of length   
                       the total number of constraints. 
                     - ``results['normxiJ']`` : :math:`L^\infty` norms of the nullspace step :math:`\xi_J`  
                       (without normalization)
                     - ``results['normxiJ_save']`` : :math:`L^\infty` norms of the nullspace step :math:`\xi_J` 
                       after every renormalization. 
                     - ``results['eps']`` : threshold values estimated for the constraints and used to take into    
                       account constraints which as close to be saturated or violated (these are those smaller than   
                       these ``eps`` values).
                     - ``results['tolerance']``: estimation of an uncertainty bound on the
                       constraints under which these can expect to be
                       satisfied. The uncertainty bound is computed thanks to the formula:
                       :math:`\texttt{tolerance} = ||\texttt{DC}||_1 \texttt{params['dt']}` 
                       where ``DC`` is the Jacobian matrix of the constraints.
                     - ``results['s']``: the optimization path length
                       ``[s(x0),s(x1),...,s(xn)]``
                       with :math:`\texttt{s(x(t))}=\int_0^t ||\texttt{x}'(\tau)||_2 d\tau`.
                        
                     It is possible to save only a limited number of previous iterations, see        
                     :ref:`params['save_only_N_iterations'] <paramsN_iterations>`.



    **Algorithm parameters**

    - ``params['alphaJ']``: (default 1) scaling coefficient for the null space
      step :math:`\xi_J` decreasing the objective function. :math:`\alpha_J(t)` is computed     
      such that :math:`||\alpha_J(t)\xi_J(x(t))||_{\infty}\leqslant \texttt{params['alphaJ']}`.

    - ``params['alphaC']``   : (default 1) scaling coefficient for the Gauss Newton
      direction :math:`\xi_C` decreasing the violation of the constraints.  
      :math:`\alpha_C(t)` is computed     
      such that :math:`||\alpha_C(t)\xi_C(x(t))||_{\infty}\leqslant \texttt{params['alphaC']}`.

    - ``params['alphas']``: (default ``None``) a list ``[alpha1,...,alpha_{p+q}]``  
      or a function   
      ``lambda x : alphas(x)`` where
      ``alphas(x)=[alpha1,...,alpha_{p+q}]`` is a list of weights 
      scaling the Gauss Newton direction :math:`\xi_C` for
      each of the constraints at the point ``x``.   
      Here ``p`` is the number of equality constraints and ``q`` the number of inequality 
      constraints.
      This parameter can be  a function to allow
      for cases where the number of constraints may depend on ``x`` (for instance 
      if ``x`` is a discretization of an infinite dimensional object) or if it 
      is important to change the weights depending on ``x``.

    - ``params['debug']``: Tune the level of verbosity of the output (default 0).   
      The higher this parameter and the more verbose is the output.
      Use ``param['debug']=-1`` to display only the final result
      Use ``param['debug']=-2`` to remove any output.

    - ``params['dt']``: (default : `0.1`). Time-step
      used for the discretization of :eq:`nlspace_ode` with the Forward Euler method.     

    - ``params['itnormalisation']``: (default 1)    
      the 
      null space step :math:`\xi_J` is normalized until the iteration
      ``params['itnormalisation']``.  
      More precisely, :math:`\alpha_J(t)` is computed such that 
      :math:`||\alpha_J(t)\xi_J(x(t))||_{\infty}=\texttt{params['alphaJ']}`
      until the iteration number ``params['itnormalisation']`` (after    
      this equality becomes an inequality).
        
      Increasing this parameter allows to go faster during the first iteration but may also results in a    
      two large time step.

    - ``params['K']``: tunes the distance at which inactive inequality constraints
      are "felt"; this avoids oscillations around active constraints. 
      nstraints are "felt" from a distance ``K*params['dt']`` from the boundary.


    - ``params['maxit']``: Maximal number of iterations (default : 4000)

    - ``params['maxtrials']``: (default 3) number of trials in between time steps
      until the finite difference check is acceptable (see ``params['tol_finite_diff']``) 

    - ``params['normalisation_norm']`` : the norm used to normalize the
      descent direction (default is :math:`L^\infty` encoded by ``numpy.inf``).

    - ``params['dual_norm']`` : the dual norm of the norm provided by
      ``params['normalisation_norm']``.

    - ``params['normalize_tol']`` : (default: -1) if >= 0 ,
      then :math:`\xi_J` is normalized such that    :math:`||\alpha_J(t)\xi_J(x(t))||_{\infty}=\texttt{params['alphaJ']}`
      every time the set of active constraints changes
      (in addition to the first ``params['itnormalisation']`` normalization iterations).    
      This allows to speed up the optimizer when there are a small number of constraints.
      The
      value of this parameter can be set to a strictly positive tolerance
      (e.g. 1e-7) to normalize only when a substantial discontinuity occurs
      in the multipliers. If set to a negative value, then no normalization
      is performed when the active set changes.
        
      .. note::

         Turning this parameter on is not recommended if there are many constraints as     
         too many normalizations 
         can make the convergence to the local optimum unstable.   

    - ``params['tol']``: (default 1e-5). The algorithm stops when
      :math:`||x_{n+1}-x_n||_{2}<\texttt{params['tol']}\times \texttt{params['dt']}``
      or after ``params['maxit']`` iterations.

    - ``params['tol_finite_diff']`` : (default 0.15) a warning message will be  
      displayed if ``params[`debug`]=1`` and if the finite difference checks

      .. math ::    
        
         \begin{aligned}  
         |J(x_{n+1})-J(x_n) - \texttt{DJ(dx)}| & \leqslant \texttt{params['tol\_finite\_diff']}\times ||\texttt{dx}||_2 \\ 

         | C(x_{n+1})-C(x_n) - dC(dx) |  & \leqslant \texttt{params['tol\_finite\_diff']}\times ||\texttt{dx}||_2
         \end{aligned}
        
      fails, where :math:`C` is the vector of constraints. Useful to check the  
      correct implementation of the sensitivities.

    - ``params['tol_qp']`` : (default 1e-10) the tolerance for the qp solver.

    - ``params['tol_cg']`` : (default 1e-15) the tolerance for scipy iterative linear solvers   
      (e.g. `lsqr <https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.lsqr.html>`_).

    - ``params['show_progress_qp']`` : (default ``False``) If set to ``True``,  
      the output of
      qp solver will be displayed between iterations.

    - ``params['qp_solver']`` : ``osqp`` (default) ``cvxopt`` or ``qpalm``. 
      The quadratic programming solver used to
      solve the dual problem,   
      either `OSQP <https://osqp.org/docs/examples/setup-and-solve.html>`_, 
      `CVXOPT <https://cvxopt.org/>`__ or `QPALM <https://github.com/kul-optec/QPALM>`_.

    - ``params['qp_solver_options']`` : a dictionary of options that can be passed  
      to the quadratic programming solver. Check    
      the documentation of  
      `OSQP <https://osqp.org/docs/interfaces/solver_settings.html>`__ and
      `CVXOPT <https://cvxopt.org/userguide/coneprog.html>`__.

      .. _paramsN_iterations:

    - ``params['save_only_N_iterations']`` : (default ``None``) if set to a integer ``N>0``,    
      then the routine will only save the last ``N`` iterates     
      in ``results['x']``.  The past iterates are set to ``None``, for saving memory.
      By default, the routine save all iterates, which can be costly for design variables     
      containing a lot of information such as large arrays.

    - ``params['save_only_Q_constraints']`` : (default ``None``)    
      this parameter is taken into account only   
      if ``params['save_only_N_iterations']`` is set to a positive integer. Then,    
      if set to an integer ``Q``,     
      the routine will save data related to only the first ``Q`` equality and ``Q``   
      inequality constraints past the last ``N`` iterations.  
      This enables to save memory when the number of constraints is large, e.g.   
      in the case of bound constraints.
        
    - ``params['start']`` : (default ``None``)  
      If set to an integer, and if a dictionary of previous ``results`` is provided,    
      the optimizer will restart at iteration ``params['start']``. 

    - ``params['CFL']`` : TODO
        
    - ``params['method_xiC']`` : (default ``linear_system``). Specify the method    
        for computing ``xiC``. By default, will try to solve the linear system  
        :math:`DC \xi_C=C`. If set to ``qp``, then null space will solve a quadratic  
        program instead (may be more robust against degenerate constraints). 

    """

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

        A = problem.inner_product(x)

        dJ = problem.dJ(x)
        dG = problem.dG(x)
        dH = problem.dH(x)

        J, G, H, dJ, dG, dH, C, dC, n, p, q = pack_constraints(J, G, H, dJ, dG, dH)

        prevmuls = muls.copy()

        # Null space direction xiJ and range space direction xiC
        tic()
        xiJ, xiC, eps, tildeEps, muls = get_xiJ_xiC(J, G, H, dJ, dG, dH, A,  
                                               h=params['K']*params['dt'],        
                                               alphas = params['alphas'](x),
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
        dx = -AJ*xiJ-AC*xiC
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
    return abstract_results.implementation()
