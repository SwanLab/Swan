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
    params = {'alphaC': 1, 'debug': 0, 'alphaJ': 1, 'dt': 0.1, 'maxtrials': 1}
    params.update(other_params)
    results = nlspace_solve(basicProblem(), params)
    return results

def solve_basic_problem_equalized(**other_params):
    params = {'alphaC': 1, 'debug': 0, 'alphaJ': 1, 'dt': 0.1, 'maxtrials': 1}
    params.update(other_params)
    resultsEqualized = nlspace_solve(
        EqualizedOptimizable(basicProblem()), params)
    resultsEqualized['x'] = [list(x[0])+list(x[1])
                             for x in resultsEqualized['x']]
    return resultsEqualized

def restart_basic_problem(results, **other_params):
    resultsCopy = results.copy()
    for key in results:
        resultsCopy[key] = results[key].copy()[:100]
    params = {'alphaC': 1, 'debug': 0, 'alphaJ': 1, 'dt': 0.1, 'maxtrials': 1}
    params.update(other_params)
    return nlspace_solve(basicProblem(), params, resultsCopy)
        
def main(**options):
    import nullspace_optimizer.examples.basic_examples.utils.draw as draw

    results = solve_basic_problem(**options)
    print(f"Optimum : {results['x'][-1]}")
    resultsEqualized = solve_basic_problem_equalized(**options)
    print(f"Method of slack ended in {len(resultsEqualized['J'])} iterations.")
    print(f"Nullspace method ended in {len(results['J'])} iterations.")

    input("\nWill test restarting, press any key")
    resultsNew = restart_basic_problem(results, **options)

    print("\nResults after restart:")
    for key in resultsNew:
        print("{0:<10} before restart: \t".format(key), len(results[key]),
              " \t after restart:\t", len(resultsNew[key]))

    print("")
    print("Optimum :")
    print(results['x'][-1])
    print("Comparison with restarting : ")
    print(resultsNew['x'][-1])

    draw.ion()
    draw.drawProblem(basicProblem(), XLIM=[0.2, 2.8], YLIM=[
                     0.2, 2.8], resolution=200)
    draw.drawData(resultsEqualized, 'EQUALIZED', 'green', x0=True)
    draw.drawData(results, 'NLSPACE', 'blue')

    draw.figure()
    draw.drawMuls(results, 'NLSPACE')
    draw.legend()

    draw.figure()
    draw.drawJ(results)
    draw.drawJ(resultsEqualized, 'EQUALIZED', linestyle='--')
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
    if '--qpalm' in sys.argv:
        options.update(qp_solver='qpalm')
    if '--qpverbose' in sys.argv:
        options.update(show_progress_qp=True)
    if '--K' in sys.argv:   
        options.update(K=0)
    main(**options)
