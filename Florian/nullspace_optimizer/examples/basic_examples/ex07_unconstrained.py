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


class unconstrainedProblem(EuclideanOptimizable):
    def x0(self):
        return [1.5, 2.25]

    def J(self, x):
        return x[1]**2+0.1*x[0]**2

    def dJ(self, x):
        return [0.2*x[0], 2*x[1]]


def run_problems(**other_params):
    params = {'alphaC': 1, 'debug': 0, 'alphaJ': 1, 'dt': 1, 'maxtrials': 1}
    params.update(other_params)
    return nlspace_solve(unconstrainedProblem(), params)

def main():
    import nullspace_optimizer.examples.basic_examples.utils.draw as draw

    results = run_problems()

    print(f"Nullspace method ended in {len(results['J'])} iterations.")

    print("")
    print("Optimum :")
    print(results['x'][-1])

    draw.ion()
    draw.drawProblem(unconstrainedProblem(), XLIM=[-3, 3], YLIM=[-3, 3],
                resolution=200)
    draw.drawData(results, 'NLSPACE', 'blue')

    draw.figure()
    draw.drawJ(results)
    draw.legend()

    draw.show()

    input("\nPress any key to close all plots")
    draw.close('all')

if __name__ == "__main__":
    main()
