# Copyright 2018-2019 CNRS, Ecole Polytechnique and Safran.
#
# This file is part of nullspace_optimizer.
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

from nullspace_optimizer import EuclideanOptimizable
from nullspace_optimizer import nlspace_solve


class LinearDependentOptimizable(EuclideanOptimizable):
    def x0(self):
        return [0.5, 2.]

    def J(self, x):
        return (x[0]**2+x[1]**2)*0.5

    def H(self, x):
        return [-x[0],
                -2*x[0],
                -x[1]+1]

    def dJ(self, x):
        return x

    def dH(self, x):
        return [[-1,0],
                [-2,0],
                [0,-1]]

linear_pb = LinearDependentOptimizable()

def run_problem():  
    params = {'dt': 0.1, 'debug': -1}
    results = nlspace_solve(linear_pb, params)
    return results

def main():
    import nullspace_optimizer.examples.basic_examples.utils.draw as draw
    results = run_problem()
    draw.ion()
    draw.drawProblem(linear_pb, [-0.5, 1.5], [-0.5, 2.5])
    draw.drawData(results, label='x', color = f'C{1}', x0=True, xfinal=True, initlabel=None)
    input("Press any key to close the figures")


if __name__ == "__main__":
    main()

