import numpy as np
import numpy.testing as npt
    
from nullspace_optimizer.examples.advanced_optimizable_examples import ex00_bound_constraints,\
    ex02_filter,\
    ex03_symbolic_svanberg,\
    ex04_symbolic_svanberg,\
    ex05_basic_example_symbolic        

def test_ex00_bound_constraints():   
    results = ex00_bound_constraints.run(N=10, plot=False)
    assert results['J'][-1]<=-42    
    assert np.max(results['x'][-1])<=1
    assert np.min(results['x'][-1])>=0

def test_ex01_bound_constraints_sparse():   
    results = ex00_bound_constraints.run(N=100, plot=False,  
                                         maxit=10, qp_solver='osqp')
    assert results['J'][-1]<=-447
    assert np.max(results['x'][-1])<=1
    assert np.min(results['x'][-1])>=0

def test_ex01_bound_constraints_cvxopt():   
    results = ex00_bound_constraints.run(N=100, plot=False, qp_solver='cvxopt',  
                                         maxit=10)
    assert results['J'][-1]<=-447   
    assert np.max(results['x'][-1])<=1
    assert np.min(results['x'][-1])>=0
        

def test_ex01_save_only_N_iterations():
    results = ex00_bound_constraints.run(N=100, plot=False,  
                                         maxit=10, qp_solver='osqp')
    for x in results['x'][:-1]: 
        assert x is None    
    for key in ['H','eps','muls','tolerance']:  
        for data in results[key][:-1]:  
            assert len(data) == 5
        
def test_ex02_filter():
    results = ex02_filter.solve_basic_problem()
    actual = ex02_filter.filter(results['x'][-1])
    expected = np.asarray([1.825589497529696, 0.5477682695635483])
    npt.assert_allclose(actual, expected, rtol=1.0e-5)
        
def test_ex03_symbolic():
    results = ex03_symbolic_svanberg.solve_P1()
    expected = np.array([5.        , 6.01601754, 5.30917567, 4.4943317 , 3.50147811,
                         2.1526566 ]) 
    npt.assert_allclose(results['x'][-1], expected)

def test_ex04_symbolic():
    results = ex04_symbolic_svanberg.solve_P1()
    expected = np.array([0.        , 1.41163046, 0.3770739 ]) 
    npt.assert_allclose(results['x'][-1], expected)

def test_ex05_basicpb_symbolic():
    results = ex05_basic_example_symbolic.solve_basic_problem()
    expected = np.asarray([1.82570174, 0.54773459])
    npt.assert_allclose(results['x'][-1], expected)
    
def test_cvxopt_osqp():
    case = ex00_bound_constraints.Problem(10)   
    x = case.x0()   
    J = case.J(x)
    G = case.G(x)
    H = case.H(x)   
    dJ = case.dJ(x)
    dG = case.dG(x)     
    dH = case.dH(x) 
    A = case.inner_product(x)
    from nullspace_optimizer.optimizers.nullspace.utils import get_tilde, pack_constraints  
    J, G, H, dJ, dG, dH, C, dC, n, p, q = pack_constraints(J, G, H, dJ, dG, dH)
    tilde = get_tilde(C, p)
    from nullspace_optimizer.optimizers.nullspace.osqp_interface import OsQpSolver
    from nullspace_optimizer.optimizers.nullspace.cvx_interface import CvxQpSolver
    qpSolverOsQp = OsQpSolver()
    qpSolverCVX = CvxQpSolver()
    resOsQp = qpSolverOsQp(A,dJ,dC[tilde,:],p)
    resCVX = qpSolverCVX(A,dJ,dC[tilde,:],p)
    objOsqp=0.5*resOsQp[:n].T @ A @ resOsQp[:n]+dJ @ resOsQp[:n]
    objCVX=0.5*resCVX[:n].T @ A @ resCVX[:n]+dJ @ resCVX[:n]
    assert np.isclose(objOsqp,objCVX,rtol=1e-5)
    constraintCVX = A @ resCVX[:n]-dC[tilde,:].T @ resCVX[n:]
    constraintOsQp = A @ resOsQp[:n]-dC[tilde,:].T @ resOsQp[n:]
    assert np.allclose(constraintCVX,0,atol=1e-15)
    assert np.allclose(constraintOsQp,0,atol=1e-15)
    mulsCVX = np.zeros(p+q)
    mulsOsqp = np.zeros(p+q)
    mulsCVX[tilde] = resCVX[n:]
    mulsOsqp[tilde] = resOsQp[n:]
    assert(min(mulsCVX)>-1e-15)
    assert(min(mulsOsqp)>-1e-15)
        
def test_nullspace():
    case = ex00_bound_constraints.Problem(100)   
    from nullspace_optimizer.optimizers.nullspace.utils import get_xiJ_xiC, pack_constraints,\
                                                               get_eps
    x = case.x0()   
    J = case.G(x)
    G = case.G(x)
    H = case.H(x)   
    dJ = case.dJ(x)
    dG = case.dG(x)     
    dH = case.dH(x) 
    A = case.inner_product(x)
    J, G, H, dJ, dG, dH, C, dC, n, p, q = pack_constraints(J, G, H, dJ, dG, dH)
    (eps, tildeEps) = get_eps(C, dC, p, h=0.1)
    xiJ, _, _, _, _ = get_xiJ_xiC(J, G, H, dJ, dG, dH, A, h=0.1)
    if p>0:
        assert np.isclose(dC[:p,:] @ xiJ,0,atol=1e-15) 
    if dC[tildeEps,:][p:,:].size > 0:
        assert np.min(dC[tildeEps,:][p:,:]@xiJ)>-1e-15
            
            
