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





path = "NonLinearNull/Summum/"

exports = FreeFemRunner(path+"MeshMBB.edp").execute()
Th = exports['Th']

## POSTPROCESS
data  = np.load(path+"Results/MBBTOrigPlotting.npz",allow_pickle=True)
iterT = data['it']
CompT = data['c']
VT  = data['v']
thick = data['thick']

data  = np.load(path+"Results/MBBT5Plotting.npz",allow_pickle=True)
iterT5 = data['it']
CompT5 = data['c']
VT5  = data['v']
thick5 = data['thick']

data = np.load(path+"Results/MBBTOrigPlotting.npz")
fig0, ax0 = plt.subplots()
x = data["xF"]
u = P1Function(Th,x<=0)
u.plot(fig=fig0,ax=ax0)
plt.pause(0.05)

fig, axes = plt.subplots(1, 3, figsize=(10, 4))

axes[0].plot(iterT, CompT, color='b', label='1 geo it')
axes[0].plot(iterT5, CompT5, color='r', label='5 geo it')
axes[0].legend()
axes[0].set_xlabel('Iter')
axes[0].set_ylabel('Compliance')
axes[0].grid(True, linestyle='--', alpha=0.6)

axes[1].plot(iterT, VT, color='b', label='1 geo it')
axes[1].plot(iterT5, VT5, color='r', label='5 geo it')
axes[1].legend()
axes[1].set_xlabel('Iter')
axes[1].set_ylabel('Volume constraint')
axes[1].grid(True, linestyle='--', alpha=0.6)

axes[2].plot(iterT, thick, color='b', label='1 geo it')
axes[2].plot(iterT5, thick5, color='r', label='5 geo it')
axes[2].legend()
axes[2].set_xlabel('Iter')
axes[2].set_ylabel('Max. Thick. constraint')
axes[2].grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.show()
plt.pause(0.1)

a = 1