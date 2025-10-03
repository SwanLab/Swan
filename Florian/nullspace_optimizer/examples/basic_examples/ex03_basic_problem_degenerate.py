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

class basicProblemDegenerate(EuclideanOptimizable):
    def x0(self):
        return [1.5, 2.25]

    def J(self, x):
        return x[1]+0.3*x[0]

    def dJ(self, x):
        return [0.3, 1]

    def H(self, x):
        return [-x[1]+1.0/x[0], -(3-x[0]-x[1]), 0*x[0]]

    def dH(self, x):
        return [[-1.0/x[0]**2, -1], [1, 1], [0,0]]


def run_basic_problem(**other_params):
    params = {'alphaC': 1, 'debug': 0, 'alphaJ': 1, 'dt': 0.1, 'maxtrials': 1}
    params.update(other_params)
    results = nlspace_solve(basicProblemDegenerate(), params)
    return results


def run_basic_problem_equalized(**other_params):
    params = {'alphaC': 1, 'debug': 0, 'alphaJ': 1, 'dt': 0.1, 'maxtrials': 1}
    params.update(other_params)
    resultsEqualized = nlspace_solve(
        EqualizedOptimizable(basicProblemDegenerate()), params)
    resultsEqualized['x'] = [list(x[0])+list(x[1])
                             for x in resultsEqualized['x']]
    return resultsEqualized


def main(**options):
    results = run_basic_problem(**options)
    resultsEqualized = run_basic_problem_equalized(**options)
    print(f"Method of slack ended in {len(resultsEqualized['J'])} iterations.")
    print(f"Nullspace method ended in {len(results['J'])} iterations.")

    print("")
    print("Optimum :")
    print(results['x'][-1])
    print("Optimum with Equalized ")
    print(resultsEqualized['x'][-1][:2])

if __name__ == "__main__":
    import sys
    options = dict()
    if '--osqp' in sys.argv:
        options.update(qp_solver='osqp')
    if '--qpverbose' in sys.argv:
        options.update(show_progress_qp=True)
    main(**options)
