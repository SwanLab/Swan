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

# This test case is an adaptation of topopt_cholmod.py  from
# A 200 LINE TOPOLOGY OPTIMIZATION CODE BY NIELS AAGE
# AND VILLADS EGEDE JOHANSEN, JANUARY 2013
from nullspace_optimizer import EuclideanOptimizable, bound_constraints_optimizable,\
                                memoize, filtered_optimizable
from nullspace_optimizer import nlspace_solve
import numpy as np
from nullspace_optimizer.examples.basic_examples.utils import draw
import cvxopt
import cvxopt.cholmod
import scipy.sparse as sp
import matplotlib.pyplot as plt
from matplotlib import colors
import sys
#element stiffness matrix
def lk():
    E=1
    nu=0.3
    k=np.array([1/2-nu/6,1/8+nu/8,-1/4-nu/12,-1/8+3*nu/8,-1/4+nu/12,-1/8-nu/8,nu/6,1/8-3*nu/8])
    KE = E/(1-nu**2)*np.array([     
        [k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
        [k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
        [k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
        [k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
        [k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
        [k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
        [k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
        [k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]] ]);
    return (KE)

def deleterowcol(A, delrow, delcol):
    # Assumes that matrix is in symmetric csc form !
    m = A.shape[0]
    keep = np.delete (np.arange(0, m), delrow)
    A = A[keep, :]
    keep = np.delete (np.arange(0, m), delcol)
    A = A[:, keep]
    return A    

# Max and min stiffness
Emin=1e-9
Emax=1.0
    
def init(**kwargs):
    global nelx, nely, volfrac, rmin, penal, H, Hs, ndof, KE, iK, jK, fixed, free
    global edofMat, f, dofs
    # Default input parameters
    nelx=kwargs.get('nelx',180)
    nely=kwargs.get('nely',60)
    volfrac=kwargs.get('volfrac',0.4)
    rmin=kwargs.get('rmin',5.4)
    penal=kwargs.get('penal',3.0)

    if __name__=='__main__':
        if len(sys.argv)>1: nelx   =int(sys.argv[1])
        if len(sys.argv)>2: nely   =int(sys.argv[2])
        if len(sys.argv)>3: volfrac=float(sys.argv[3])
        if len(sys.argv)>4: rmin   =float(sys.argv[4])
        if len(sys.argv)>5: penal  =float(sys.argv[5])

    ndof = 2*(nelx+1)*(nely+1)
    dofs = np.arange(ndof).reshape((nelx+1,nely+1,2))
    edofMat = np.zeros((nelx,nely,8),dtype=int) 
    edofMat[:,:,0] = dofs[:-1,1:,0]
    edofMat[:,:,1] = dofs[:-1,1:,1]
    edofMat[:,:,2] = dofs[1:,1:,0]
    edofMat[:,:,3] = dofs[1:,1:,1]
    edofMat[:,:,4] = dofs[1:,:-1,0]
    edofMat[:,:,5] = dofs[1:,:-1,1]
    edofMat[:,:,6] = dofs[:-1,:-1,0]
    edofMat[:,:,7] = dofs[:-1,:-1,1]
    edofMat=edofMat.reshape((nelx*nely,8))
    # FE: Build the index vectors for the for coo matrix format.
    KE=lk()

    # Construct the index pointers for the coo format
    iK = np.kron(edofMat,np.ones((8,1))).flatten()
    jK = np.kron(edofMat,np.ones((1,8))).flatten()    

    # BC's and support
    fixed=np.union1d(dofs[0,:,0],dofs[-1,-1,1])
    free=np.setdiff1d(dofs.reshape(ndof),fixed)

    
    # Filter: Build (and assemble) the index+data vectors for the coo matrix format
    nfilter=int(nelx*nely*((2*(np.ceil(rmin)-1)+1)**2))
    iH = np.zeros(nfilter)
    jH = np.zeros(nfilter)
    sH = np.zeros(nfilter)
    cc=0
    for i in range(nelx):
        for j in range(nely):
            row=i*nely+j
            kk1=int(np.maximum(i-(np.ceil(rmin)-1),0))
            kk2=int(np.minimum(i+np.ceil(rmin),nelx))
            ll1=int(np.maximum(j-(np.ceil(rmin)-1),0))
            ll2=int(np.minimum(j+np.ceil(rmin),nely))
            for k in range(kk1,kk2):
                for l in range(ll1,ll2):
                    col=k*nely+l
                    fac=rmin-np.sqrt(((i-k)*(i-k)+(j-l)*(j-l)))
                    iH[cc]=row
                    jH[cc]=col
                    sH[cc]=np.maximum(0.0,fac)
                    cc=cc+1
    # Finalize assembly and convert to csc format
    H=sp.coo_matrix((sH,(iH,jH)),shape=(nelx*nely,nelx*nely)).tocsc()    
    Hs=sp.diags(1/np.asarray(H.sum(1)).flatten())
    

@memoize()
def filter(x):  
    return Hs @ (H @ x) 
    
@memoize()
def diff_filter(x,v):
    return (H @ (Hs @ v.T)).T

@memoize()
def solve_state(x, f = None):
    dc = np.ones(nely*nelx)
    ce = np.ones(nely*nelx)

    # Setup and solve FE problem
    sK=((KE.flatten()[np.newaxis]).T*(Emin+x**penal*(Emax-Emin))).flatten(order='F')
    K = sp.coo_matrix((sK,(iK,jK)),shape=(ndof,ndof)).tocsc()
    # Remove constrained dofs from matrix and convert to coo
    K = deleterowcol(K,fixed,fixed).tocoo()
    # Solve system 
    K = cvxopt.spmatrix(K.data,K.row.astype(int),K.col.astype(int))
        
    if f is None:   
        f=np.zeros((ndof,1))
        f[1,0]=-1
    B = cvxopt.matrix(f[free,:])
    cvxopt.cholmod.linsolve(K,B)
    u=np.zeros((ndof,f.shape[1]))
    u[free,:]=np.array(B)[:,:] 

    ce = np.einsum('ija,jk,ika->ia',u[edofMat,:],KE,u[edofMat,:])
    obj = (Emin + (x**penal*(Emax-Emin))).dot(ce)
    dc = np.multiply((-penal*x**(penal-1)*(Emax-Emin))[:,np.newaxis],ce).T

    return (obj,dc)
    
@bound_constraints_optimizable(l=0,u=1)
@filtered_optimizable(filter, diff_filter)
class TO_problem(EuclideanOptimizable):
    def __init__(self, plot=None):
        self.plot = plot
        if self.plot:
            # Initialize plot and plot the initial design
            plt.ion() # Ensure that redrawing is possible
            self.fig,self.ax = plt.subplots()

    def x0(self):
        return volfrac * np.ones(nely*nelx,dtype=float)

    def J(self, x):
        (obj,dc)=solve_state(x)
        return obj[0]
        
    def dJ(self, x):
        (obj,dc)=solve_state(x)
        return dc[0,:]

    def G(self,x):
        return [np.sum(x)/(nelx*nely)-volfrac]

    def dG(self,x):
        dv = np.ones(nelx*nely)/(nelx*nely)
        return dv

    def accept(self, params, results):
        if self.plot:
            x = results['x'][-1]

            if not hasattr(self,'im'):
                self.im = self.ax.imshow(-x.reshape((nelx,nely)).T, cmap='gray',\
                                    interpolation='none',norm=colors.Normalize(vmin=-1,vmax=0))
                self.ax.axis('off')
                self.fig.show()
            else:
                self.im.set_array(-x.reshape((nelx,nely)).T)
                self.fig.canvas.draw()
            plt.pause(0.01)
                

            
optimization_params = {"dt": 0.3, 'alphaJ': 1,
                       'alphaC': 1, 
                       'maxtrials': 1,
                       'maxtrials':3,
                       'show_progress_qp':0,
                       'normalize_tol':-1,
                       'itnormalisation': 30,
                       'tol_finite_diff':10,  
                       'save_only_N_iterations':1,    
                       'save_only_Q_constraints':5,     
                       'qp_solver':'qpalm',
                       'tol_qp':1e-9,       
                       'maxit':150}
                
if __name__ == "__main__":
    init()
    case=TO_problem(plot=True)
    results = nlspace_solve(case,
                            optimization_params)
    plt.figure()
    draw.drawMuls(results)
    plt.figure()
    draw.drawJ(results)
    input("Press any key")
