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

from nullspace_optimizer import nlspace_solve
import matplotlib.pyplot as plt
from nullspace_optimizer.examples.basic_examples.utils import draw

from nullspace_optimizer.examples.advanced_optimizable_examples.ex00_bound_constraints import run

def main(N=100,**params):
    plt.ion()
    results = run(N,**params)
    plt.figure()
    draw.drawMuls(results)
    plt.figure()
    draw.drawJ(results)

    input("Press any key")

if __name__=="__main__":
    main(100, qp_solver='osqp', plot=True, debug=10)
