from nullspace_optimizer import nlspace_solve
from nullspace_optimizer.examples.topopt_examples import ex00_compliance as ex00    
    
    
    
if __name__=="__main__":
    ex00.init(nelx=150, nely = 50, volfrac=0.01, rmin=1.9, penal = 3)
    case = ex00.TO_problem(plot = True) 
    results = nlspace_solve(case, {"dt":0.02, "itnormalisation":30, "maxit":150,    
                                   'qp_solver':'qpalm','tol_qp':1e-8})

