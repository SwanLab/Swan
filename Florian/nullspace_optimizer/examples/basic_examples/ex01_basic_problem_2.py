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

class basicProblem2(EuclideanOptimizable):
    def x0(self):
        return [1.5, 2.25]

    def J(self, x):
        return (x[0]-2)**2+(x[1]-2)**2

    def dJ(self, x):
        return [2*(x[0]-2), 2*(x[1]-2)]

    def H(self, x):
        return [-x[1]+1.0/x[0], -(3-x[0]-x[1])]

    def dH(self, x):
        return [[-1.0/x[0]**2, -1], [1, 1]]
        
def solve_basic_problem_2(**other_params):
    params = {'alphaC': 0.2, 'alphaJ': 1, 'dt': 0.1}
    params.update(other_params)
    return nlspace_solve(basicProblem2(), params)

def solve_basic_problem_2_equalized(**other_params):
    params = {'alphaC': 0.2, 'alphaJ': 1, 'dt': 0.1}
    params.update(other_params)
    resultsEqualized = nlspace_solve(
        EqualizedOptimizable(basicProblem2()), params)
    resultsEqualized['x'] = [list(x[0])+list(x[1])
                             for x in resultsEqualized['x']]
    return resultsEqualized

def main(**options):
    import nullspace_optimizer.examples.basic_examples.utils.draw as draw
    results = solve_basic_problem_2(**options)
    resultsEqualized = solve_basic_problem_2_equalized(**options)


    print("")
    print("Optimum :")
    print(results['x'][-1])
    print(f"Method of slack ended in {len(resultsEqualized['J'])} iterations.")
    print(f"Nullspace method ended in {len(results['J'])} iterations.")

    draw.ion()
    draw.drawProblem(basicProblem2(), XLIM=[0.04, 3], YLIM=[
                     0.2, 2.7], resolution=200)
    draw.drawData(resultsEqualized, 'EQUALIZED',
                  'green', loc='lower left', x0=True)
    draw.drawData(results, 'NLSPACE', 'blue', loc='lower left')

    draw.figure()
    draw.drawMuls(results, 'NLSPACE')
    draw.legend()

    draw.figure()
    draw.drawJ(results)
    draw.drawJ(resultsEqualized, 'EQUALIZED', linestyle='--')
    draw.legend()

    draw.figure()
    draw.drawC(results)
    draw.drawC(resultsEqualized, 'EQUALIZED', linestyle='--')
    draw.legend()
        
    draw.show() 

    input("Press any key to close all plots")

if __name__ == "__main__":
    import sys
    options = dict()
    if '--debug' in sys.argv:
        options.update(debug=5)
    if '--osqp' in sys.argv:
        options.update(qp_solver='osqp')
    if '--qpverbose' in sys.argv:
        options.update(show_progress_qp=True)
    main(**options)
