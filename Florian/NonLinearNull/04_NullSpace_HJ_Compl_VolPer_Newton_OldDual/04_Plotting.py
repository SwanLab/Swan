## PACKAGES
import sys
sys.path.append("/home/joseantonio/Documentos/GitHub/Swan/Florian"), \

from nullspace_optimizer import EuclideanOptimizable,\
bound_constraints_optimizable, memoize, filtered_optimizable
import numpy as np
from nullspace_optimizer.optimizable import Optimizable
from nullspace_optimizer.inout import tic, toc   
from nullspace_optimizer import inout as io   
from nullspace_optimizer.optimizers.nullspace.utils import compute_norm, get_xiJ_xiC, get_tilde, pack_constraints
from nullspace_optimizer.optimizers.utils import OptimizationResults, check_params
from nullspace_optimizer.optimizers.nullspace import utils
import cvxopt
import cvxopt.cholmod
import scipy.sparse as sp
import matplotlib.pyplot as plt
from matplotlib import colors
from pyfreefem import FreeFemRunner
from pymedit import P1Function





path = "NonLinearNull/04_NullSpace_HJ_Compl_VolPer_Newton_OldDual/"

exports = FreeFemRunner(path+"04_Mesh.edp").execute()
Th = exports['Th']

## POSTPROCESS
data  = np.load(path+"04_ResultCaseOriginal.npz")
iter1 = data['it']
Comp1 = data['c']
Vol1  = data['v']
Per1  = data['p']

data  = np.load(path+"04_ResultCaseNewton.npz")
iter2 = data['it']
Comp2 = data['c']
Vol2  = data['v']
Per2  = data['p']

data = np.load(path+"04_ResultCaseNewton.npz")
fig0, ax0 = plt.subplots()
x = data["xF"]
u = P1Function(Th,x<=0)
u.plot(fig=fig0,ax=ax0)
plt.pause(0.05)

fig, axes = plt.subplots(1, 3, figsize=(10, 4))

axes[0].plot(iter1, Comp1, color='b', label='Original')
axes[0].plot(iter2, Comp2, color='r', label='Almost Newton')
axes[0].legend()
axes[0].set_xlabel('Iter')
axes[0].set_ylabel('Compliance')
axes[0].grid(True, linestyle='--', alpha=0.6)

axes[1].plot(iter1, Vol1, color='b', label='Original')
axes[1].plot(iter2, Vol2, color='r', label='Almost Newton')
axes[1].legend()
axes[1].set_xlabel('Iter')
axes[1].set_ylabel('Volume constraint')
axes[1].grid(True, linestyle='--', alpha=0.6)

axes[2].plot(iter1, Per1, color='b', label='Original')
axes[2].plot(iter2, Per2, color='r', label='Almost Newton')
axes[2].legend()
axes[2].set_xlabel('Iter')
axes[2].set_ylabel('Perimeter constraint')
axes[2].grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.show()
plt.pause(0.1)

a = 1