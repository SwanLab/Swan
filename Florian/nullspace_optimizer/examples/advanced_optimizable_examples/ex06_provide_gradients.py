import numpy as np
from nullspace_optimizer import EuclideanOptimizable, nlspace_solve, memoize
    
theta = np.pi/4
R = np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
A = R @ np.array([[10,0],[0,1]]) @ R.T
Ainv = np.linalg.inv(A)

class basicProblem(EuclideanOptimizable):
    def x0(self):
        return [1.5, 2.25]

    def J(self, x):
        return x[1]+0.3*x[0]

    @memoize(func_name="dJ",debug=2)
    def dJ(self, x):
        return [0.3, 1]

    def H(self, x):
        return [-x[1]+1.0/x[0], -(3-x[0]-x[1])]

    @memoize(func_name="dH",debug=2)
    def dH(self, x):
        return np.array([[-1.0/x[0]**2, -1], [1, 1]])

    def dJT(self, x):
        return Ainv @ self.dJ(x)

    def dHT(self, x):
        return Ainv @ (self.dH(x).T)   

def solve_basic_problem(**other_params):
    params = {'dt': 0.3,    
              'maxtrials': 1, 'provide_gradients': True, 'debug':3}
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
