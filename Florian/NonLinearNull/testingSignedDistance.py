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
import skfmm
from matplotlib.collections import PolyCollection
from matplotlib import colors

import matlab.engine
eng = matlab.engine.start_matlab()
eng.addpath(eng.genpath('C:/Users/JOSE A. TORRES/Documents/GitHub/Swan'), nargout=0)
sig = eng.SignedDistance(nargout=1)
 
def to_matlab(x):
    import matlab
    return matlab.double(np.atleast_2d(x).tolist())

def from_matlabVector(x):
    import numpy as np
    return np.fromiter((v[0] for v in x),dtype=np.float64)


[phi,nx,ny,dx] = eng.getPhi(sig,nargout=4)
phi = from_matlabVector(phi)
nx = np.int64(nx)
ny = np.int64(ny)
phi = np.transpose(np.reshape(phi,(ny,nx)))
dx = float(np.double(dx))
U = skfmm.distance(phi,dx)
map = to_matlab(U)

eng.plotSignedDistance(sig,map,nargout=0)

end=True