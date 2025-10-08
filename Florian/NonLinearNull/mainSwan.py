import sys
sys.path.append("C:/Users/JOSE A. TORRES/Documents/GitHub/Swan/Florian"), \

from nullspace_optimizer import EuclideanOptimizable,\
bound_constraints_optimizable, memoize, filtered_optimizable,\
nlspace_solve
import numpy as np
import cvxopt
import cvxopt.cholmod
import scipy.sparse as sp
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib import colors
import matlab.engine
eng = matlab.engine.start_matlab()
eng.addpath(eng.genpath('C:/Users/JOSE A. TORRES/Documents/GitHub/Swan'), nargout=0)
tutorial = eng.TutorialToPythonDensitySetting(nargout=1)
 
def init(**kwargs):
# Initialized global variables
    global nelx, nely, volfrac, rmin, penal
    global H, Hs, ndof, KE, iK, jK, fixed, free
    global edofMat, f, dofs
 
 
# Default input parameters
    nelx = kwargs.get("nelx", 180)
    nely = kwargs.get("nely", 60)
    volfrac = kwargs.get("volfrac", 0.4)
    rmin = kwargs.get("rmin", 5.4)
    penal = kwargs.get("penal", 3.0)

def to_matlab(x):
    import matlab
    return matlab.double(np.atleast_2d(x).tolist())

def from_matlab(x):
    import numpy as np
    return np.fromiter((v[0] for v in x),dtype=np.float64) #np.array(x)


# Definition of the optimization problem
@bound_constraints_optimizable(l=0, u=1)
class TO_problem(EuclideanOptimizable):
    def x0(self):
        xI = eng.getInitialGuess(tutorial)
        return from_matlab(xI)

    def J(self, x):
        xM = to_matlab(x)
        [j,dj] = eng.computeCost(tutorial,xM,nargout=2)
        obj = np.float64(j)
        return obj

    def dJ(self, x):
        xM = to_matlab(x)
        [j,dj] = eng.computeCost(tutorial,xM,nargout=2)
        dc = from_matlab(dj)
        return dc

    def G(self, x):
        xM = to_matlab(x)
        [j,dj] = eng.computeConstraint(tutorial,xM,nargout=2)
        obj = np.array([np.float64(j)],dtype=np.float64)
        return obj
    def dG(self, x):
        xM = to_matlab(x)
        [j,dj] = eng.computeConstraint(tutorial,xM,nargout=2)
        dc = from_matlab(dj)
        return dc

    def accept(self, params, results):
        # Plot the design at every iteration
        [xM,patch] = eng.filterRho(tutorial,nargout=2)
        x = from_matlab(xM)
        if not hasattr(self, "fig"):
            self.fig, self.ax = plt.subplots()
            patches = np.array(patch)
            self.collection = PolyCollection(
                patches,
                edgecolors='none',
                facecolors='gray',
                alpha=1.0
            )
            self.ax.add_collection(self.collection)
            self.ax.set_aspect('equal')
            self.ax.axis('off')
            self.fig.show()
            self.ax.set_xlim(0, 2)
            self.ax.set_ylim(0, 1)
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
        cmap = plt.get_cmap('gray')
        facecolors = cmap(1- x)
        self.collection.set_facecolor(facecolors)
        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()
        plt.pause(0.01)

# Optimization parameters
optimization_params = {"dt": 0.05,
                        "itnormalisation": 50,
                        "save_only_N_iterations": 1,
                        "save_only_Q_constraints": 5,
                        "alphaJ": 1.25,
                        "alphaC": 1,
                        "maxit": 150}
# Initialize and solve the TO problem
init()
case = TO_problem()
results = nlspace_solve(case, optimization_params) 

it = results['it']
J = results['J']
G = results['G']

fig, axes = plt.subplots(1, 2, figsize=(10, 4))

axes[0].plot(it, J, color='b')
axes[0].set_xlabel('Iter')
axes[0].set_ylabel('Compliance')
axes[0].grid(True, linestyle='--', alpha=0.6)

axes[1].plot(it, G, color='b')
axes[1].set_xlabel('Iter')
axes[1].set_ylabel('Volume constraint')
axes[1].grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.show()