import sys
sys.path.append("/home/joseantonio/Documentos/GitHub/Swan/Florian"), \
    # C:/Users/JOSE A. TORRES/Documents/GitHub/Swan/Florian

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
 
# element stiffness matrix
def lk():
    E = 1
    nu = 0.3
    k = np.array([1/2-nu/6, 1/8+nu/8,-1/4-nu/12,-1/8+3*nu/8,
        -1/4+nu/12,-1/8-nu/8, nu/6, 1/8-3*nu/8])
    KE = E/(1-nu**2)*np.array([
        [k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
        [k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
        [k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
        [k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
        [k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
        [k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
        [k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
        [k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]]])
    return (KE)
 
def deleterowcol(A, delrow, delcol):
    # Delete row and columns of a sparse csc matrix A
    m = A.shape[0]
    keep = np.delete(np.arange(0, m), delrow)
    A = A[keep, :]
    keep = np.delete(np.arange(0, m), delcol)
    A = A[:, keep]
    return A
 
# Max and min stiffness
Emin = 1e-9
Emax = 1.0
 
# Assemble FEM matrices and filter
def init(**kwargs):
# Initialized global variables
    global nelx, nely, volfrac, rmin, penal
    global HFilt, Hs, ndof, KE, iK, jK, fixed, free
    global edofMat, f, dofs
 
 
# Default input parameters
    nelx = kwargs.get("nelx", 180)
    nely = kwargs.get("nely", 60)
    volfrac = kwargs.get("volfrac", 0.4)
    rmin = kwargs.get("rmin", 5.4)
    penal = kwargs.get("penal", 3.0)

# Degrees of Freedom matrix: edofMat[i,j,:] contains the
# DOFs associated to the square element (i,j)
    ndof = 2*(nelx+1)*(nely+1)
    dofs = np.arange(ndof).reshape((nelx+1, nely+1, 2))
    edofMat = np.zeros((nelx, nely, 8), dtype=int)
    edofMat[:, :, 0] = dofs[:-1, 1:, 0]
    edofMat[:, :, 1] = dofs[:-1, 1:, 1]
    edofMat[:, :, 2] = dofs[1:, 1:, 0]
    edofMat[:, :, 3] = dofs[1:, 1:, 1]
    edofMat[:, :, 4] = dofs[1:, :-1, 0]
    edofMat[:, :, 5] = dofs[1:, :-1, 1]
    edofMat[:, :, 6] = dofs[:-1, :-1, 0]
    edofMat[:, :, 7] = dofs[:-1, :-1, 1]
    edofMat = edofMat.reshape((nelx*nely, 8))

# FE: Build the index vectors for the for coo matrix format.
    KE = lk()
# Construct the index pointers for the coo format
    iK = np.kron(edofMat, np.ones((8, 1))).flatten()
    jK = np.kron(edofMat, np.ones((1, 8))).flatten()

# BC’s and support
    fixed = np.union1d(dofs[0, :, 0], dofs[-1,-1, 1])
    free = np.setdiff1d(dofs.reshape(ndof), fixed)

# Filter: Build (and assemble) the index+data vectors
# for the coo matrix format
    nfilter = int(nelx*nely*((2*(np.ceil(rmin)-1)+1)**2))
    iH = np.zeros(nfilter)
    jH = np.zeros(nfilter)
    sH = np.zeros(nfilter)
    cc = 0
    for i in range(nelx):
        for j in range(nely):
            row = i*nely+j
            kk1 = int(np.maximum(i-(np.ceil(rmin)-1), 0))
            kk2 = int(np.minimum(i+np.ceil(rmin), nelx))
            ll1 = int(np.maximum(j-(np.ceil(rmin)-1), 0))
            ll2 = int(np.minimum(j+np.ceil(rmin), nely))
            for k in range(kk1, kk2):
                for l in range(ll1, ll2):
                    col = k*nely+l
                    fac = rmin-np.sqrt(((i-k)*(i-k)+(j-l)*(j-l)))

                    iH[cc] = row
                    jH[cc] = col
                    sH[cc] = np.maximum(0.0, fac)
                    cc = cc+1
# Finalize assembly and convert to csc format
    HFilt = sp.coo_matrix((sH, (iH, jH)),
        shape=(nelx*nely, nelx*nely)).tocsc()
    Hs = sp.diags(1/np.asarray(HFilt.sum(1)).flatten())

# f is the unit load vector at the bottom right corner
    f = np.zeros((ndof, 1))
    f[1, 0] =-1

# Compute compliance and its sensitivity
@memoize()
def solve_state(x):
    dc = np.ones(nely*nelx)
    ce = np.ones(nely*nelx)

    # Setup and solve FE problem
    sK = ((KE.flatten()[np.newaxis]).T*(Emin+x**penal*(Emax-Emin)))\
    .flatten(order="F")
    K = sp.coo_matrix((sK, (iK, jK)), shape=(ndof, ndof)).tocsc()
    # Remove constrained dofs from matrix and convert to coo
    K = deleterowcol(K, fixed, fixed).tocoo()
    # Solve system
    K = cvxopt.spmatrix(K.data, K.row.astype(int), K.col.astype(int))
    B = cvxopt.matrix(f[free, :])
    cvxopt.cholmod.linsolve(K, B)
    u = np.zeros((ndof, f.shape[1]))
    u[free, :] = np.array(B)[:, :]

    # compliance
    ce = np.einsum("ija,jk,ika->ia", u[edofMat, :], KE, u[edofMat, :])
    obj = (Emin + (x**penal*(Emax-Emin))).dot(ce)

    # sensitivity of compliance with respect to x
    dc = np.multiply((-penal*x**(penal-1)*(Emax-Emin))[:, np.newaxis],
    ce).T

    # If not multi-objective
    if len(obj)==1:
        return (obj[0], dc[0,:])
    # else
    return (obj, dc)

# Density filter
@memoize()
def filter(x):
    return Hs @ (HFilt @ x)


# Sensitivity of the filter
@memoize()
def diff_filter(x, v):
    return (HFilt @ (Hs @ v.T)).T

# Definition of the optimization problem
@bound_constraints_optimizable(l=0, u=1)
@filtered_optimizable(filter, diff_filter)
class TO_problem(EuclideanOptimizable):
    def x0(self):
        return volfrac * np.ones(nely*nelx, dtype=float)

    def J(self, x):
        (obj, dc) = solve_state(x)
        return obj

    def dJ(self, x):
        (obj, dc) = solve_state(x)
        return dc

    def G(self, x):
        return [np.sum(x)/(nelx*nely)-volfrac]

    def dG(self, x):
        dv = np.ones(nelx*nely)/(nelx*nely)
        return dv

    def accept(self, params, results):
        # Plot the design at every iteration
        x = results["x"][-1]
        if not hasattr(self, "im"):
            plt.ion() # Ensure that redrawing is possible
            self.fig, self.ax = plt.subplots()
            self.im = self.ax.imshow(-x.reshape((nelx, nely)).T,\
            cmap="gray",
            interpolation="none",
            norm=colors.Normalize(vmin=-1, vmax=0))
            self.ax.axis("off")
            self.fig.show()
        else:
            self.im.set_array(-x.reshape((nelx, nely)).T)
            self.fig.canvas.draw()
        plt.pause(0.01)

# Optimization parameters
params = {"dt": 0.05,
          "itnormalisation": 50,
          "save_only_N_iterations": 1,
          "save_only_Q_constraints": 5,
          "alphaJ": 1.25,
          "alphaC": 1,
          "maxit": 10}
# Initialize and solve the TO problem
init()
problem:Optimizable = TO_problem()







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



results = abstract_results.implementation()
iter = results['it']
Per  = results['J']
Vol  = results['G']

fig, axes = plt.subplots(1, 2, figsize=(10, 4))

axes[0].plot(iter, Per, color='b')
axes[0].set_xlabel('Iter')
axes[0].set_ylabel('Perimeter')
axes[0].grid(True, linestyle='--', alpha=0.6)

axes[1].plot(iter, Vol, color='b')
axes[1].set_xlabel('Iter')
axes[1].set_ylabel('Volume constraint')
axes[1].grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.show()

a = 1