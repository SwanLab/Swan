from nullspace_optimizer import EqualizedOptimizable, nlspace_solve
from nullspace_optimizer.examples.topopt_examples.ex00_compliance import TO_problem, init
from nullspace_optimizer.examples.basic_examples.utils import draw
import matplotlib.pyplot as plt
    
if __name__=="__main__":
    init()
    case = EqualizedOptimizable(TO_problem(plot=__name__=="__main__"))
    results = nlspace_solve(case,
                            {"dt": 0.2, 'alphaJ': 1,
                             'alphaC': 1, 
                             'maxtrials': 1,
                             'maxtrials':3,
                             'itnormalisation': 30,
                             'tol_finite_diff':10,  
                             'save_only_N_iterations':1,    
                             'save_only_Q_constraints':5,
                             'maxit':200})
    plt.figure()
    draw.drawMuls(results)
    plt.figure()
    draw.drawJ(results)
    input("Press any key")


