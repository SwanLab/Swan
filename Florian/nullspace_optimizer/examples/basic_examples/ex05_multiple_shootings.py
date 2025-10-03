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

from nullspace_optimizer import EuclideanOptimizable, nlspace_solve

class problemeSimple(EuclideanOptimizable):
    def __init__(self, xinit):
        self.xinit = xinit

    def x0(self):
        return self.xinit

    def J(self, x):
        return (x[0]+1)**2+(x[1]+1)**2

    def G(self, x):
        return []

    def H(self, x):
        return [x[0]**2+x[1]**2-1**2,
                x[0]+x[1]-1,
                -x[1]-0.7]

    def dJ(self, x):
        return [2*(x[0]+1), 2*(x[1]+1)]

    def dG(self, x):
        return []

    def dH(self, x):
        return [[2*x[0], 2*x[1]],
                [1, 1],
                [0, -1]]

def run_problems(**other_params):
    xinits = ([1.25, 0], [0, -1.2], [1, -1], [0.7, 1.2], [-1, 1])
    problems = [problemeSimple(xinit=x0) for x0 in xinits]
    params = {'dt': 0.1, 'alphaC': 1, 'alphaJ': 1,   
              'maxtrials': 1, 'debug': -1}
    params.update(other_params)
    return [nlspace_solve(pb, params) for pb in problems]

def main():
    import nullspace_optimizer.examples.basic_examples.utils.draw as draw

    results = run_problems()

    draw.ion()
    draw.drawProblem(problemeSimple([0, 0]), [-1.5, 1.5], [-1.5, 1.5])
    for i, r in enumerate(results):
        draw.drawData(r, f'x{i+1}', f'C{i}', x0=True, xfinal=True, initlabel=None)

    draw.figure()
    draw.drawMuls(results[0], 'x1')
    draw.legend()
    draw.title('Lagrange multipliers for trajectory x1')
    draw.show()
    draw.figure()
    draw.drawMuls(results[1], 'x2')
    draw.legend()
    draw.title('Lagrange multipliers for trajectory x2')
    draw.show()
    input("Press any key to close all plots")
    draw.close('all')

if __name__ == "__main__":
    main()
