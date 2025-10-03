from nullspace_optimizer.optimizers.MMA import mma_lib, mma_solve
from nullspace_optimizer.optimizers.OC import oc_solve
with_ipopt = True
try:
    from nullspace_optimizer.optimizers.IPOPT import ipopt_solve, ipopt_lower_inf, ipopt_upper_inf
except:
    print("IPOPT seems not installed. Containing without.")
    with_ipopt = False
from nullspace_optimizer.examples.topopt_examples import ex00_compliance
from nullspace_optimizer import EqualizedOptimizable
ex00_compliance.init(nelx=90, nely=30) #Low resolution
    

from nullspace_optimizer.examples.topopt_examples import ex01_compliance_equalized,\
    ex03_compliance_MMA,\
    ex07_heat_nlspace,\
    ex08_heat_MMA,  \
    ex09_heat_OC, \
    ex10_multiple_load,\
    ex11_multiple_load_MMA
    
if with_ipopt:  
    from nullspace_optimizer.examples.topopt_examples import ex04_compliance_IPOPT, \
       ex12_multiple_load_IPOPT 

ex07_heat_nlspace.init(30)
    

from nullspace_optimizer import nlspace_solve
import numpy as np


def test_ex00_compliance(): 
    case = ex00_compliance.TO_problem()
    results = nlspace_solve(case,
                            {"dt": 0.2, 'alphaJ': 1,
                             'alphaC': 1, 
                             'maxtrials': 1,
                             'maxtrials':3,
                             'itnormalisation': 30,
                             'tol_finite_diff':10,  
                             'save_only_N_iterations':1,    
                             'save_only_Q_constraints':5,   
                             'qp_solver':'osqp',
                             'maxit':5})
    assert results['J'][-1] <= 676

def test_ex00_compliance_cvxopt(): 
    case = ex00_compliance.TO_problem()
    results = nlspace_solve(case,
                            {"dt": 0.2, 'alphaJ': 1,
                             'alphaC': 1, 
                             'maxtrials': 1,
                             'maxtrials':3,
                             'itnormalisation': 30,
                             'tol_finite_diff':10,  
                             'save_only_N_iterations':1,    
                             'save_only_Q_constraints':5,   
                             'qp_solver':'cvxopt',
                             'maxit':5})
    assert results['J'][-1] <= 760
        
def test_ex01_compliance_equalized():
    case = EqualizedOptimizable(ex00_compliance.TO_problem())   
    results_osqp = nlspace_solve(case,
                                 {"dt": 0.2, 'alphaJ': 1,
                                  'alphaC': 1, 
                                  'maxtrials': 1,
                                  'maxtrials':3,
                                  'itnormalisation': 30,
                                  'tol_finite_diff':10,  
                                  'save_only_N_iterations':1,    
                                  'save_only_Q_constraints':5,
                                  'maxit':5, 
                                  'qp_solver':'osqp'})
    results_cvxopt = nlspace_solve(case,
                                   {"dt": 0.2, 'alphaJ': 1,
                                    'alphaC': 1, 
                                    'maxtrials': 1,
                                    'maxtrials':3,
                                    'itnormalisation': 30,
                                    'tol_finite_diff':10,  
                                    'save_only_N_iterations':1,    
                                    'save_only_Q_constraints':5,
                                    'maxit':5, 
                                    'qp_solver':'cvxopt'})
    assert results_osqp['J'][-1] <= 760
    assert results_osqp['J'][-1] == results_cvxopt['J'][-1]
        
def test_ex02_oc(): 
    case = ex00_compliance.TO_problem()
    results = oc_solve(case, 0, 1,
                          {'move':0.2,  
                           'maxit':5})
    assert results['J'][-1]<=760
        
def test_ex03_mma():
    case = ex03_compliance_MMA.TO_problem()
    results = ex03_compliance_MMA.mma_solve(case, 0, 1,     
                                            {'move': 0.2,    
                                             'maxit':5,
                                             'a0':1,
                                             'c':10000*np.ones((1,1))})
    assert results['J'][-1]<=730

if with_ipopt:
    def test_ex04_ipopt():
        case = ex04_compliance_IPOPT.TO_problem()
        results = ipopt_solve(case, 0, 1,     
                              {'maxit':5})
        assert results['J'][-1]<=765
        
def test_ex05_tiny_volfrac():   
    ex00_compliance.volfrac = 0.01 #Low resolution
    ex00_compliance.rmin = 1.9 #Low resolution
    case = ex00_compliance.TO_problem()
    results = nlspace_solve(case, {"dt":0.01,"itnormalisation":30,      
                                   "qp_solver":"qpalm", 
                                   "tol_qp": 1e-6,
                                   "maxit":5})
    assert results['J'][-1]<=3.86e7 and abs(results['G'][-1][0])<1e-4
    
def test_ex06_tiny_volfrac_MMA():   
    ex00_compliance.volfrac = 0.01 #Low resolution
    ex00_compliance.rmin = 1.9 #Low resolution
    case = ex03_compliance_MMA.TO_problem() 
    results = mma_solve(case, 0, 1, {'move':0.2,  
                                     'fix':True,
                                     'c':1000,
                                     'maxit':5,   
                                     'tight_move':True})
    assert results['J'][-1]<=3.98e7 and results['H'][-1][0]<0
    
def test_ex07_heat_nlspace():    
    case = ex07_heat_nlspace.Heat_TO()
    params = dict(dt=0.3, itnormalisation=30, maxit=5)  
    results = nlspace_solve(case, params)
    assert results['J'][-1] <= 771 
    
def test_ex08_heat_MMA():    
    case = ex08_heat_MMA.Heat_TO()  
    params = dict(move=0.3, maxit=5)
    results = mma_solve(case, 0, 1, params)
    assert results['J'][-1]<=853
    
def test_ex09_heat_OC(): 
    case = ex09_heat_OC.Heat_TO()  
    params = dict(move=0.3,maxit=5)
    results = oc_solve(case, 0, 1, params)
    assert results['J'][-1]<=850
    
def test_ex10_multiple_load():  
    ex00_compliance.init(nelx=90,nely=30)
    ex10_multiple_load.init()
    case = ex10_multiple_load.TO_problem()
    params = dict(dt=0.2, maxit=5, itnormalisation=10)
    results = nlspace_solve(case, params)
    assert results['J'][-1] <= 7

def test_ex11_multiple_load_MMA():  
    ex00_compliance.init(nelx=90,nely=30)
    ex10_multiple_load.init()
    case = ex11_multiple_load_MMA.TO_problem()
    params = dict(move=0.1, maxit=5)
    l = [0]*ex00_compliance.nelx*ex00_compliance.nely + [-1e5]    
    u = [1]*ex00_compliance.nelx*ex00_compliance.nely + [1e5]
    results = mma_solve(case, l, u, params)
    assert results['J'][-1] <= 5

if with_ipopt:
    def test_ex12_multiple_load_IPOPT():  
        ex00_compliance.init(nelx=90,nely=30)
        ex10_multiple_load.init()
        case = ex12_multiple_load_IPOPT.TO_problem()
        params = dict(maxit=5)
        l = [0]*ex00_compliance.nelx*ex00_compliance.nely + [ipopt_lower_inf]    
        u = [1]*ex00_compliance.nelx*ex00_compliance.nely + [ipopt_upper_inf]
        results = ipopt_solve(case, l, u, params)
        assert results['J'][-1] <= 9
    
