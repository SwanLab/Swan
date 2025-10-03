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

import matplotlib.pyplot as plt
import matplotlib as mp
import numpy as np
from nullspace_optimizer import EuclideanOptimizable

mp.rcParams['contour.negative_linestyle'] = 'solid'
mp.rcParams['lines.linewidth'] = 2
mp.rcParams['axes.labelsize'] = 16
mp.rcParams['xtick.labelsize'] = 16
mp.rcParams['ytick.labelsize'] = 16
mp.rcParams['legend.fontsize'] = 15
mp.rcParams['figure.autolayout'] = True

with_plot = True

def drawProblem(problem: EuclideanOptimizable, XLIM, YLIM, resolution=100):
    """Draw the admissible domain of a 2-d optimization problem.

    Inputs
    ------

    problem :   an instance of EuclideanOptimizable in two dimensions
                (problem.n=2)
    XLIM, YLIM: dimensions of the box in which the problem will be drawn:
                XLIM=[xmin,xmax], YLIM=[ymin,ymax]
    resolution: (default 100) number of grid points in each directions.
    """
    x = np.linspace(XLIM[0], XLIM[1], resolution)
    y = np.linspace(YLIM[0], YLIM[1], resolution)
    X, Y = np.meshgrid(x, y)
    ZJ = problem.J([X, Y])
    ax = plt.gca()
    ax.patch.set_alpha(0.3)
    CS = plt.contour(X, Y, ZJ, alpha=0.5, colors='k', linewidths=1)
    plt.clabel(CS, inline=1, inline_spacing=2, fontsize=10)
    H = problem.H([X,Y])
    if H:
        H = np.maximum.reduce(H)
        plt.contourf(X, Y, H, levels=[2*np.min(H), 0.0, 2*np.max(H)],
                     colors=['w', 'silver', 'silver'], alpha=0.7)
        plt.contour(X, Y, H, levels=[0.0], colors=['r'])
    ax.set_aspect('equal')
    plt.xlabel(r'$x_0$', fontsize=20)
    plt.ylabel(r'$x_1$', fontsize=20)


def drawData(results, label, color, loc='upper right', x0=None, linewidth=2,
             xfinal=True, **kwargs):
    """Draw optimization path provided by the nlspace_solve function

    Inputs
    ------

    results  : output of a nlspace_solve run
    label    : label of optimization path
    color    : color of optimization path and final iterate
    loc      : parameter for the legend
    x0       : (default None) if True, then the initialization is shown
    linewidth: line width of optimization path
    xfinal   : if True, then a marker is plotted at the final point
    """
    x1 = [x[0] for x in results['x']]
    x2 = [x[1] for x in results['x']]
    if x0:
        plt.plot(x1[0], x2[0], '.', color=color, markersize=10,
                 label=kwargs.pop('initlabel', 'Initialization'))
    if xfinal:
        plt.plot(x1[-1], x2[-1], 'X', color=color, markersize=10,
                 linewidth=3, alpha=0.2)
    plt.plot(x1, x2, linestyle=kwargs.pop('linestyle','--'), markersize=kwargs.pop('markersize',5), 
             label=label,
             color=color, linewidth=linewidth, **kwargs)
    ax = plt.gca()
    ax.legend(loc=loc, fontsize=10)
    ax.set_aspect('equal')
    ax.autoscale(tight=True)


def drawMuls(results, title='Muls', linewidth=2, path_length=False, **kwargs):
    """Draw Lagrange multipliers provided by the nlspace_solve function
    """
    muls = list(zip(*results['muls']))
    maxit = kwargs.pop('maxit', len(results['it']))
    if path_length:
        abscissa = results['s'][:maxit]
    else:
        abscissa = list(range(maxit))
    for i, mul in enumerate(muls):
        plt.plot(abscissa, mul[:maxit], linewidth=linewidth,
                 label=r'$\mu_'+str(i)+r'$'+f' - {title}', **kwargs)
        if i>10:
            break
    if path_length:
        plt.xlabel('s', fontsize=16)
    plt.title(title)
    plt.legend()
    plt.gca().legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2)


def drawJ(results, label='NLSPACE', linewidth=2, path_length=False, **kwargs):
    """Draw Objective function provided by the nlspace_solve function"""
    maxit = kwargs.pop('maxit', len(results['it']))
    if path_length:
        abscissa = results['s'][:maxit]
    else:
        abscissa = list(range(maxit))
    plt.plot(abscissa, results['J'][:maxit],
             color='C0', label='J - '+label, linewidth=linewidth, **kwargs)
    if path_length:
        plt.xlabel('s', fontsize=16)
    plt.gca().legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2)


def drawC(results, label='NLSPACE', linewidth=2, path_length=False, **kwargs):
    """Draw Constraint functions provided by the nlspace_solve function"""
    maxit = kwargs.pop('maxit', len(results['s']))
    if path_length:
        abscissa = results['s'][:maxit]
    else:
        abscissa = range(maxit)
    if len(results['H'][0]) > 0:
        H1 = [x[0] for x in results['H']]
        H2 = [x[1] for x in results['H']]
    else:
        H1 = [gi[0]-zi[2]**2/2 for (gi, zi)
              in zip(results['G'], results['x'])]
        H2 = [gi[1]-zi[3]**2/2 for (gi, zi)
              in zip(results['G'], results['x'])]
    plt.plot(abscissa, H1[:maxit],
             label=r'$h_1$ '+f'({label})', linewidth=linewidth, color='C0',
             **kwargs)
    plt.plot(abscissa, H2[:maxit],
             label=r'$h_2$ '+f'({label})', linewidth=linewidth, color='C1',
             **kwargs)
    if path_length:
        plt.xlabel('s', fontsize=16)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2)

def ion():
    plt.ion()

def figure(*args, **kwargs):
    plt.figure(*args, **kwargs)

def title(*args, **kwargs):
    plt.title(*args, **kwargs)

def legend(*args, **kwargs):
    plt.legend(*args, **kwargs)

def remove_legend():
    plt.gca().get_legend().remove()

def show(*args, **kwargs):
    plt.show(*args, **kwargs)

def close(*args, **kwargs):
    plt.close(*args, **kwargs)
