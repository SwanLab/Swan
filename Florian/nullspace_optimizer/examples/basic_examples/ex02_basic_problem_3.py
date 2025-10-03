# This file is part of nullspace_optimizer.
#   
# This file has been modified from a version    
# released under 
# Copyright 2018-2019 CNRS, Ecole Polytechnique and Safran.
#
# nullspace_optimizer is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# nullspace_optimizer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# A copy of the GNU General Public License is included below.
# For further information, see <http://www.gnu.org/licenses/>.

from nullspace_optimizer import nlspace_solve, EqualizedOptimizable, EuclideanOptimizable

class parabProblem(EuclideanOptimizable):
    def x0(self):
        return [3, 3]

    def J(self, x):
        return x[1]**2+(x[0]+3)**2

    def dJ(self, x):
        return [2*(x[0]+3), 2*x[1]]

    def H(self, x):
        return [-x[0]**2+x[1], -x[1]-x[0]-2]

    def dH(self, x):
        return [[-2*x[0], 1], [-1, -1]]


def solve_problem_parab(**other_params):
    params = {'alphaC': 0.2, 'alphaJ': 1,  'dt': 0.4, 'qp_solver':'cvxopt'}
    params.update(other_params)
    return nlspace_solve(parabProblem(), params)


def solve_problem_parab_equalized(**other_params):
    params = {'alphaC': 0.2, 'alphaJ': 1,  'dt': 0.4}
    params.update(other_params)
    resultsEqualized = nlspace_solve(
        EqualizedOptimizable(parabProblem()), params)
    resultsEqualized['x'] = [list(x[0])+list(x[1])
                             for x in resultsEqualized['x']]
    return resultsEqualized

def solve_problem_parab_osqp(**other_params):
    params = {'alphaC': 0.2, 'alphaJ': 1,  'dt': 0.4, 'qp_solver': 'osqp'}
    params.update(other_params)
    return nlspace_solve(parabProblem(), params)

def main(**options):
    import nullspace_optimizer.examples.basic_examples.utils.draw as draw
    from nullspace_optimizer.inout import drawHistories 
    results = solve_problem_parab(**options)
    input("Continue")
    resultsEqualized = solve_problem_parab_equalized(**options)
    input("Continue")
    resultsOSQP = solve_problem_parab_osqp(**options)
    input("Continue")

    draw.ion()
    draw.drawProblem(parabProblem(), XLIM=[-3.5, 5], YLIM=[-1.3, 3.2])
    draw.drawData(resultsEqualized, 'EQUALIZED',
                  'green', loc='lower right', x0=True)
    draw.drawData(resultsOSQP, 'NLSPACE (OSQP)', 'orange',
                  loc='lower right', x0=True, linestyle='-')
    draw.drawData(results, 'NLSPACE', 'blue', loc='lower right')
        
    print("")
    print("Optimum :")
    print(results['x'][-1])
    print(f"Method of slack ended in {len(resultsEqualized['J'])} iterations.")
    print(
        f"Nullspace method (CVXOPT) ended in {len(results['J'])} iterations.")
    print(
        f"Nullspace method (OSQP) ended in {len(resultsOSQP['J'])} iterations.")

    draw.figure()
    draw.drawMuls(results, 'NLSPACE (CVXOPT)')
    draw.drawMuls(resultsOSQP, 'NLSPACE (OSQP)', linestyle='--')
    draw.legend()

    draw.figure()
    draw.drawJ(resultsOSQP, linestyle=':')
    draw.drawJ(results)
    draw.drawJ(resultsEqualized, 'EQUALIZED', linestyle='--')
    draw.legend()

    draw.figure()
    draw.drawC(resultsOSQP, linestyle=':')
    draw.drawC(resultsEqualized, 'EQUALIZED', linestyle='--')
    draw.legend()

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

