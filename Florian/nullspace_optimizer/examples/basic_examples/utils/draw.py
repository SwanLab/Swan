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

try:
    from nullspace_optimizer.examples.basic_examples.utils.draw_pyplot import *

except ImportError:
    from nullspace_optimizer import EuclideanOptimizable

    with_plot = False

    def drawProblem(problem: EuclideanOptimizable, XLIM, YLIM, resolution=100):
        pass

    def drawMuls(results, title="Muls", linewidth=2, **kwargs):
        pass

    def drawJ(results, label='NLSPACE', linewidth=2, **kwargs):
        pass

    def drawC(results, label='NLSPACE', linewidth=2, **kwargs):
        pass

    def drawData(results, label, color, loc='upper right', x0=None, linewidth=2,
                 xfinal=True, **kwargs):
        pass

    def ion():
        pass

    def figure(*args, **kwargs):
        pass

    def title(*args, **kwargs):
        pass

    def legend(*args, **kwargs):
        pass

    def remove_legend():
        pass

    def show(*args, **kwargs):
        pass

    def close(*args, **kwargs):
        pass
