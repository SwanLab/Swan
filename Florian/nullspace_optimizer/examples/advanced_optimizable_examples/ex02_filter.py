import numpy as np
from nullspace_optimizer import filtered_optimizable, EuclideanOptimizable, nlspace_solve

def filter(x):  
    return [x[0]+x[1]**2, 2*x[1]]
        
def diff_filter(x,v):   
    return v @ np.asarray([[1,2*x[1]],[0,2]])


@filtered_optimizable(filter, diff_filter)
class basicProblem(EuclideanOptimizable):
    def x0(self):
        return [1.5, 2.25]

    def J(self, x):
        return x[1]+0.3*x[0]

    def dJ(self, x):
        return [0.3, 1]

    def H(self, x):
        return [-x[1]+1.0/x[0], -(3-x[0]-x[1])]

    def dH(self, x):
        return [[-1.0/x[0]**2, -1], [1, 1]]
        
def solve_basic_problem(**other_params):
    params = {'dt': 0.3, 'maxtrials': 1}
    params.update(other_params)
    results = nlspace_solve(basicProblem(), params)
    return results

def main(**options):
    import nullspace_optimizer.examples.basic_examples.utils.draw as draw
    results = solve_basic_problem(**options)
    print("")
    print("Optimum :")
    print(results['x'][-1])
    draw.ion()
    draw.drawProblem(basicProblem(), XLIM=[0.2, 2.8], YLIM=[
                     0.2, 2.8], resolution=200)
    draw.drawData(results, 'NLSPACE', 'blue')
    draw.show()
    input("Press any key to close all plots")
        
if __name__ == "__main__":
    import sys
    options = dict()
    if '--osqp' in sys.argv:
        options.update(qp_solver='osqp')
    if '--qpverbose' in sys.argv:
        options.update(show_progress_qp=True)
    main(**options)
