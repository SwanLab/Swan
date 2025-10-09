import numpy as np
import scipy.linalg as splinalg
import scipy.sparse as sp
import scipy.sparse.linalg as lg
from ...utils import memoize
from ...inout import tic, toc
from ... import inout as io
from ...optimizable.utils import pack_alphas
from .osqp_interface import OsQpSolver
from .cvx_interface import CvxQpSolver
from .qpalm_interface import QPALM_Solver

tol_cg = 1e-15
tol_qp = 1e-15
MAX_DIRECT_INVERSION = 100
EXACT_SOLVE = True


def compute_norm(x, norm):
    """ 
    A function for computing norm of a vector.  

    :param x: Input vector  
    :param norm: Either:    
                 * ``numpy.inf``. Then the function returns the     
                   infinity norm    
                 * ``2``. Then the function returns the averaged L^2 norm:  
                   ``sqrt(mean(x**2))``.    
                 * a custom norm function. Then the function returns ``norm(x)``."""
    if norm == np.inf:
        return np.linalg.norm(x, np.inf)
    elif norm == 2:
        return np.sqrt(np.mean(x**2))
    else:
        return norm(x)

def get_tilde(C, p, eps=0):
    tildeEps = C[p:] >= -eps
    tildeEps = np.asarray(np.concatenate(([True]*p, tildeEps)), dtype=bool)
    return tildeEps


def get_eps(C, dC, p, h, dual_norm=1):
    """ 

    :param h: distance at which to feel the constraints (typically  
              ``K*dt`` where ``dt`` is the time step and ``K=1``).
    :returns: tuple ``(eps, tildeEps)`` where ``eps`` is an array of length 
              ``dC.shape[0]`` such that at a distance ``h``, the constraint ``C[i]``    
              has value typically ``eps[i]``.
    """
    if dC.shape[0] == 0:
        return (0, np.array([], dtype=bool))

    if dual_norm == 1:
        norm1 = sp.csc_matrix.sum(abs(dC[p:, :]), 1)
    elif dual_norm == 2:
        norm1 = np.sqrt(sp.csc_matrix.multiply(
            dC[p:, :], dC[p:, :])/dC.shape[1])
    else:
        norm1 = dual_norm(dC[p:, :])

    norm1 = np.asarray(norm1).flatten()
    eps = norm1*h
    tildeEps = get_tilde(C, p, eps)
    return (eps, tildeEps)


@memoize()
def factorized(A):
    return lg.factorized(A)

@memoize()
def pack_constraints(J, G, H, dJ, dG, dH):
    G = np.asarray(G)
    H = np.asarray(H)
    C = np.concatenate((G, H))
    p = len(G)
    q = len(H)
    if p == 0:
        dG = np.empty((0, len(dJ)))
    if q == 0:
        dH = np.empty((0, len(dJ)))
    dG = sp.csc_matrix(dG)
    dH = sp.csc_matrix(dH)
    dC = sp.bmat([[dG], [dH]], format="csc")
    n = dC.shape[1]
    dJ = np.asarray(dJ)
    return J, G, H, dJ, dG, dH, C, dC, n, p, q
    
def get_qp_solver(qp_solver, qp_solver_options):
    if qp_solver == 'osqp':
        qpSolver = OsQpSolver(qp_solver_options)
    elif qp_solver == 'cvxopt':
        qpSolver = CvxQpSolver(qp_solver_options)
    elif qp_solver == 'qpalm':
        qpSolver = QPALM_Solver(qp_solver_options)
    else:
        raise Exception("Wrong qp solver provided: " +
                        qp_solver+". Available: cvxopt|osqp|qpalm")
    qpSolver.set_tolerance(tol_qp)
    return qpSolver
    
def get_xiJ_xiC(J, G, H, dJ, dG, dH, A=None, h=0.,   
                alphas = None,
                dual_norm=1,   
                qp_solver = 'osqp',     
                qp_solver_options = None,   
                method_xiC = "linear_system"):
    J, G, H, dJ, dG, dH, C, dC, n, p, q = pack_constraints(J, G, H, dJ, dG, dH)
    if A is None:
        A = sp.eye(len(dJ), format="csc")
    (eps, tildeEps) = get_eps(C, dC, p, h, dual_norm)
    qtildeEps = sum(np.where(tildeEps)[0] >= p)

    if p+qtildeEps <= MAX_DIRECT_INVERSION:
        io.display("Compute xiJ and xiC with direct inversion.", level=3)
        dJT, dCT = get_gradient_transpose(A, dJ, dC, tildeEps)
        tic()
        xiJ, muls = get_xiJ_direct(dJ, dC,dJT,dCT,p,q,qp_solver,qp_solver_options,qtildeEps,tildeEps)
        io.display("Computed null space direction in "+toc(),level=3)
        tic()
        if method_xiC == "linear_system":
            xiC = get_xiC_direct(dC, dCT, C, p, q, alphas, muls, tildeEps)
        elif method_xiC == "qp":
            xiC = get_xiC_sparse(A, C, dJ, dC, n, p, q, alphas, muls, tildeEps, method=method_xiC, qp_solver=qp_solver,
                                 qp_solver_options=qp_solver_options)
        else:   
            raise Exception("Parameter `method_xiC` should be linear_system|qp")

        io.display("Computed range space direction in "+toc(),level=3)
    else: 
        io.display("Compute xiJ and xiC with sparse mode.", level=3)
        tic()
        xiJ, muls = get_xiJ_sparse(A,dJ,dC,n,p,q,qp_solver,qp_solver_options,qtildeEps,tildeEps)
        io.display("Computed null space direction in "+toc(),level=3)
        tic()
        xiC = get_xiC_sparse(A, C, dJ, dC, n, p, q, alphas, muls, tildeEps, method=method_xiC, qp_solver=qp_solver,
                             qp_solver_options=qp_solver_options)
        io.display("Computed range space direction in "+toc(),level=3)
        
    return xiJ, xiC, eps, tildeEps, muls
    
def get_gradient_transpose(A, dJ, dC, tildeEps):    
    # Compute the gradients explicitly / do this for moderate size constraints
    solve = factorized(A)

    dJT = solve(dJ)
    dCT = np.zeros(dC.shape).T #sp.csc_matrix
    if hasattr(A, 'tocsc'):
        A = A.tocsc()
    tic()
    dJT = solve(dJ)
    for i in (x for x in range(dC.shape[0]) if tildeEps[x]):
        dCT[:, i] = solve(dC[i, :].toarray().flatten())
    io.display("Computed gradients in "+toc(), level=3)
        
    return dJT, dCT
    
def get_xiJ_direct(dJ, dC,dJT,dCT,p,q,qp_solver,qp_solver_options,qtildeEps,tildeEps):   
    muls = np.zeros(p+q)

    # Case 1 : no constraints ! Returns the gradient
    if p == 0 and qtildeEps == 0:
        return dJT, muls

    # Case 2 : no inequality constraints ! Returns the projected gradient
    if qtildeEps == 0 and p > 0:
        try:
            muls[tildeEps] = -np.linalg.solve(dC[tildeEps,:].dot(dCT[:,tildeEps]),    
                                              dC[tildeEps,:].dot(dJT))
        except:
            io.display(
                "Warning, constraints are not qualified, using "
                "pseudo-inverse.", 1, color="red")
            muls = np.zeros(p+q)
            dCdCTinv = np.linalg.pinv(dC[tildeEps, :].dot(dCT[:, tildeEps]))
            muls[tildeEps] = -dCdCTinv.dot(dC[tildeEps, :].dot(dJT))
        xiJ = dJT + dCT.dot(muls)
        return xiJ, muls
        
    # Case 3: inequality and equality constraints. Solve dual problem
    Ps = dC[tildeEps, :].dot(dCT[:, tildeEps])
    qs = dJ.dot(dCT[:, tildeEps])
    Gs = np.concatenate((np.zeros((qtildeEps, p)), np.eye(qtildeEps)), axis=1)
    l = np.zeros((qtildeEps,))
    u = np.asarray([np.inf]*qtildeEps)

    qpSolver = get_qp_solver(qp_solver,qp_solver_options)

    hat = np.asarray([True]*p+[False]*q)

    muls[tildeEps] = qpSolver.routine(Ps, qs, Gs, l, u)
    oldmuls = muls.copy()
    hat[p:] = muls[p:] > 30*tol_qp/100

    #Ignore duplicate constraints during the projection
    hat[np.logical_not(tildeEps)] = False

    # Compute null space direction xiJ
    if EXACT_SOLVE:
        try:
            muls = np.zeros(p+q)
            muls[hat] = -np.linalg.solve(dC[hat,:].dot(dCT[:,hat]),dC[hat,:].dot(dJT))
            #dCdCTinv = np.linalg.inv(dC[hat, :].dot(dCT[:, hat]))
        except Exception:
            io.display(
                "Warning, constraints are not qualified, using "
                "pseudo-inverse.", 1, color="red")
            dCdCTinv = np.linalg.pinv(dC[hat, :].dot(dCT[:, hat]))
            muls = np.zeros(p+q)
            muls[hat] = -dCdCTinv.dot(dC[hat, :].dot(dJT))

    if not np.all(muls[p:] >= -tol_qp):
        io.display("Warning, the active set has not been predicted "
                   + "correctly Using old lagrange multipliers", 1,
                   color="orange_4a")
        muls = oldmuls.copy()
    xiJ = dJT + dCT.dot(muls)
    return xiJ, muls
   
def get_xiC_direct(dC, dCT, C, p, q, alphas, muls, tildeEps):
    tilde = get_tilde(C, p)
    indicesEps = tilde
    if p+q == 0:
        hat = np.array([], dtype=bool)
    else:
        hat = np.asarray([True]*p+[False]*q)
    if not muls is None:
        hat[p:] = muls[p:] > 30*tol_qp/100

    # Compute range step direction xiC
    indicesEps[np.logical_and(tildeEps, np.logical_not(
        tilde))] = hat[np.logical_and(tildeEps, np.logical_not(tilde))]
    dCdCT = dC[indicesEps, :].dot(dCT[:, indicesEps])
        
    # Coefficients weighting each constraint
    alphas = pack_alphas(alphas, p, q)
    try:
        lamb = np.linalg.solve(dCdCT, alphas[indicesEps]*C[indicesEps])
    except Exception:
        io.display("Warning, constraints are not qualified. "
                + "Using scipy lstsq.", 1, color="red")
        lamb, _, _, _ = splinalg.lstsq(dCdCT, alphas[indicesEps]*C[indicesEps])
    xiC = dCT[:, indicesEps].dot(lamb)
    return xiC

def get_xiJ_sparse(A,dJ,dC,n,p,q,qp_solver,qp_solver_options,qtildeEps,tildeEps):
    solve = factorized(A)
    dJT = solve(dJ)

    muls = np.zeros(p+q)
    # Case 1 : no constraints ! Returns the gradient
    if qtildeEps == 0 and p == 0:
        return dJT, muls

    # Case 2 : no inequality constraints ! Returns the projected gradient
    if qtildeEps == 0 and p > 0:
        # No need to use qpSolver, compute xiJ directly
        G = sp.bmat([[A, dC[tildeEps, :].T],
                     [dC[tildeEps, :], None]],
                    format="csc")
        rhs = np.hstack(([0]*n, -dC[tildeEps, :] @ dJT))
        try:
            tic()
            res = lg.spsolve(G, rhs, permc_spec="MMD_AT_PLUS_A")
            io.display("Scipy spsolve converged in " +
                       toc(), level=5, color="magenta")
        except:
            io.display("Singular matrix, trying LSQR.",
                       level=5, color="magenta")
            tic()
            res = lg.lsqr(G, rhs, atol=tol_cg, btol=tol_cg)
            if res[1] >= 3:
                io.display(
                    "Warning: LSQR failed in the computation of the null space step", color="red", level=0)
            res = res[0]
        X = res[:n]
        muls[tildeEps] = res[n:]
        xiJ = dJT+X
        return xiJ, muls

    # Case 3: inequality and equality constraints. Solve dual problem
    qpSolver = get_qp_solver(qp_solver,qp_solver_options)
        
    res = qpSolver(A, dJ, dC[tildeEps, :], p)

    X = res[:n]
    muls[tildeEps] = res[n:]
    xiJ = dJT + X
    return xiJ, muls


def get_xiC_sparse(A, C, dJ, dC, n, p, q, alphas, muls, tildeEps, method="linear_system",    
                   qp_solver="qpalm",   
                   qp_solver_options=None):

    if A is None:
        A = sp.eye(len(dJ), format="csc")

    # Coefficients weighting each constraint
    alphas = pack_alphas(alphas, p, q)

    # First step:  identify set Istar of constraints
    #              violated or likely to saturate
    tilde = get_tilde(C, p)
    if p+q == 0:
        hat = np.array([], dtype=bool)
    else:
        hat = np.asarray([True]*p+[False]*q)
    if not muls is None:
        hat[p:] = muls[p:] > 30*tol_qp/100
            
    if method=="linear_system":
        indicesEps = tilde
        indicesEps[np.logical_and(tildeEps, np.logical_not(
            tilde))] = hat[np.logical_and(tildeEps, np.logical_not(tilde))]

        G = sp.bmat([[A, dC[indicesEps, :].T],
                     [dC[indicesEps, :], None]], format="csc")
        rhs = np.hstack(([0]*n, alphas[indicesEps]*C[indicesEps]))
        tic()
        try:
            io.display("G.shape="+str(G.shape),level=5)
            
            res = lg.spsolve(G, rhs, permc_spec="MMD_AT_PLUS_A")
            io.display("Scipy spsolve converged in " +
                       toc(), level=5, color="magenta")
            return res[:n]
        except:
            io.display("Singular matrix, trying LSQR.", level=5, color="magenta")
            tic()
            res = lg.lsqr(G, rhs, atol=tol_cg, btol=tol_cg)
            io.display("LSQR converged in "+toc(), level=5, color="magenta")
            if res[1] >= 3:
                io.display(
                    "Warning: LSQR may have failed in the computation of the "  
                    "range space step", color="red", level=0)
        return res[0][:n]
    elif method=="qp":  
        # Try new method    
        tic()
        if qp_solver == "cvxopt":   
            qpSolver = get_qp_solver("qpalm", None)
        else:
            qpSolver = get_qp_solver(qp_solver,qp_solver_options)
        P = A
        qq = np.zeros((A.shape[0],))   
        G = dC
        l = C
        u = np.zeros((p+q,))    
        u[:p]=C[:p] 
        u[p:] = np.inf
        ineq_saturated = np.asarray([False]*(p+q),dtype=bool)
        ineq_saturated[np.logical_and(tildeEps, np.logical_not(tilde))] =   \
            hat[np.logical_and(tildeEps, np.logical_not(tilde))]
        u[ineq_saturated] = C[ineq_saturated]
        xiC = qpSolver.routine(P, qq, G, l, u)
        io.display("Computed xiC using qp method (sparse mode) in "+toc(), level=3)
        return xiC
    else:   
        raise Exception("Error, parameter method should be either qp|linear_system")
        
