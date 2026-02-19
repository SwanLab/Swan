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





path = "NonLinearNull/07_NullSpace_HJ_Compl_VolMaxthick_MBB_Iterative/"

exports = FreeFemRunner(path+"07_Mesh.edp").execute()
Th = exports['Th']

## POSTPROCESS
data  = np.load(path+"07_ResultIts1Step1.npz",allow_pickle=True)
iter1 = data['it']
Comp1 = data['c']
VT1  = data['vt']
Vol1,T1 = zip(*VT1)
dmax1 = data['dmax']
muls1 = data['muls']
mulV1,mulT1 = zip(*muls1)
AJ1 = data['AJ']
AC1 = data['AC']
dH1 = data['dH']
dH1 = dH1.item()
dV1 = dH1[0]
dT1 = dH1[1]
xiC1 = data['xiC']
xiJ1 = data['xiJ']
g1 = data['g']
fOm1 = data['fOm']
fd1 = data['fd']

data  = np.load(path+"07_ResultIts5Step1.npz",allow_pickle=True)
iter5 = data['it']
Comp5 = data['c']
VT5  = data['vt']
Vol5,T5 = zip(*VT5)
dmax5 = data['dmax']
muls5 = data['muls']
mulV5,mulT5 = zip(*muls5)
AJ5 = data['AJ']
AC5 = data['AC']
dH5 = data['dH']
dH5 = dH5.item()
dV5 = dH5[0]
dT5 = dH5[1]
xiC5 = data['xiC']
xiJ5 = data['xiJ']
g5 = data['g']
fOm5 = data['fOm']
fd5 = data['fd']

data = np.load(path+"07_ResultIts1Step0.3.npz")
fig0, ax0 = plt.subplots()
x = data["xF"]
u = P1Function(Th,x<=0)
u.plot(fig=fig0,ax=ax0)
plt.pause(0.05)

fig, axes = plt.subplots(2, 4, figsize=(10, 4))

axes[0][0].plot(iter1, Comp1, color='b', label='1 geo it')
axes[0][0].plot(iter5, Comp5, color='r', label='5 geo it')
axes[0][0].legend()
axes[0][0].set_xlabel('Iter')
axes[0][0].set_ylabel('Compliance')
axes[0][0].grid(True, linestyle='--', alpha=0.6)

axes[0][1].plot(iter1, Vol1, color='b', label='1 geo it')
axes[0][1].plot(iter5, Vol5, color='r', label='5 geo it')
axes[0][1].legend()
axes[0][1].set_xlabel('Iter')
axes[0][1].set_ylabel('Volume constraint')
axes[0][1].grid(True, linestyle='--', alpha=0.6)

axes[0][2].plot(iter1, T1, color='b', label='1 geo it')
axes[0][2].plot(iter5, T5, color='r', label='5 geo it')
axes[0][2].legend()
axes[0][2].set_xlabel('Iter')
axes[0][2].set_ylabel('Max. Thick. constraint')
axes[0][2].grid(True, linestyle='--', alpha=0.6)

axes[0][3].plot(iter1, dmax1, color='b', label='1 geo it')
axes[0][3].plot(iter5, dmax5, color='r', label='5 geo it')
axes[0][3].legend()
axes[0][3].set_xlabel('Iter')
axes[0][3].set_ylabel('dmax')
axes[0][3].grid(True, linestyle='--', alpha=0.6)

axes[1][0].plot(iter1, mulV1, color='b', label='1 geo it')
axes[1][0].plot(iter5, mulV5, color='r', label='5 geo it')
axes[1][0].legend()
axes[1][0].set_xlabel('Iter')
axes[1][0].set_ylabel('Volume multiplier')
axes[1][0].grid(True, linestyle='--', alpha=0.6)

axes[1][1].plot(iter1, mulT1, color='b', label='1 geo it')
axes[1][1].plot(iter5, mulT5, color='r', label='5 geo it')
axes[1][1].legend()
axes[1][1].set_xlabel('Iter')
axes[1][1].set_ylabel('Thickness multiplier')
axes[1][1].grid(True, linestyle='--', alpha=0.6)

axes[1][2].plot(iter1, AJ1, color='b', label='1 geo it')
axes[1][2].plot(iter5, AJ5, color='r', label='5 geo it')
axes[1][2].legend()
axes[1][2].set_xlabel('Iter')
axes[1][2].set_ylabel('AJ')
axes[1][2].grid(True, linestyle='--', alpha=0.6)

axes[1][3].plot(iter1, AC1, color='b', label='1 geo it')
axes[1][3].plot(iter5, AC5, color='r', label='5 geo it')
axes[1][3].legend()
axes[1][3].set_xlabel('Iter')
axes[1][3].set_ylabel('AC')
axes[1][3].grid(True, linestyle='--', alpha=0.6)

funV1 = P1Function(Th,dV1.toarray()[0])
funT1 = P1Function(Th,dT1.toarray()[0])
funV5 = P1Function(Th,dV5.toarray()[0])
funT5 = P1Function(Th,dT5.toarray()[0])

plt.tight_layout()
plt.show()
plt.pause(0.1)

a = 1