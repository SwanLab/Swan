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





path = "NonLinearNull/09_NullSpace_HJ_SimpleGeom_VolPer_Newton_OldDual/"

exports = FreeFemRunner(path+"09_Mesh.edp").execute()
Th = exports['Th']
alpha = exports['alpha']
beta = exports['beta']

## POSTPROCESS
data  = np.load(path+"09_ResultCaseOriginal.npz")
iter1 = data['it']
Cost1 = data['c']
Vol1  = data['v']
Per1  = data['p']

data  = np.load(path+"09_ResultCaseNewton.npz")
iter2 = data['it']
Cost2 = data['c']
Vol2  = data['v']
Per2  = data['p']

data = np.load(path+"09_ResultCaseNewton.npz")
fig0, ax0 = plt.subplots()
x = data["xF"]
u = P1Function(Th,x<=0)
u.plot(fig=fig0,ax=ax0)
plt.pause(0.05)

fig, axes = plt.subplots(1, 3, figsize=(10, 4))

axes[0].plot(iter1, Cost1, color='b', label='Original')
axes[0].plot(iter2, Cost2, color='r', label='Almost Newton')
axes[0].legend()
axes[0].set_xlabel('Iter')
axes[0].set_ylabel('Geom fun')
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

runner = FreeFemRunner(path+"09_BoundaryRefinement.edp")
runner.import_variables(Th=Th,phiVal=x,alpha=alpha,lsLab=10,rInner=3)
exports = runner.execute()

Th2 = exports['Th2']
nx = exports['nx[]']
ny = exports['ny[]']

g = data["g"]
AJ = data["AJ"]
AC = data["AC"]
muls = data["muls"]
wVal = data["wVal"]
dOmVal = data["dOmVal"]
dxdOmVal = data["dxdOmVal"]
dydOmVal = data["dydOmVal"]
nullStep = data["nullStep"]
runner = FreeFemRunner(path+"09_NewtonPlotting.edp")
runner.import_variables(Th=Th,Th2=Th2,nxVal=nx,
                        nyVal=ny,g1st=g,beta=beta,lsLab=10,
                        rInner=3,rOuter=2,aJ=AJ.item(),aC=AC.item(),lam1=muls[0],lam2=muls[1],
                        wVal=wVal,dOmVal=dOmVal,
                        dxdOmVal=dxdOmVal,dydOmVal=dydOmVal)
g0 = runner.execute()['g[]']

runner = FreeFemRunner(path+"09_NewtonPlotting.edp")
runner.import_variables(Th=Th,Th2=Th2,nxVal=nx,
                        nyVal=ny,g1st=nullStep,beta=beta,lsLab=10,
                        rInner=3,rOuter=2,aJ=AJ.item(),aC=AC.item(),lam1=muls[0],lam2=muls[1],
                        wVal=wVal,dOmVal=dOmVal,
                        dxdOmVal=dxdOmVal,dydOmVal=dydOmVal)
null0 = runner.execute()['g[]']

a = 1