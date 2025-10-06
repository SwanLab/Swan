from nullspace_optimizer import EuclideanOptimizable,\
bound_constraints_optimizable, memoize, filtered_optimizable,\
nlspace_solve
import numpy as np
import cvxopt
import cvxopt.cholmod
import scipy.sparse as sp
import matplotlib.pyplot as plt
from matplotlib import colors
 
# element stiffness matrix
def lk():
    E = 1
    nu = 0.3
    k = np.array([1/2-nu/6, 1/8+nu/8,-1/4-nu/12,-1/8+3*nu/8,
        -1/4+nu/12,-1/8-nu/8, nu/6, 1/8-3*nu/8])
    KE = E/(1-nu**2)*np.array([
        [k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
        [k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
        [k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
        [k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
        [k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
        [k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
        [k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
        [k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]]])
    return (KE)
 
def deleterowcol(A, delrow, delcol):
    # Delete row and columns of a sparse csc matrix A
    m = A.shape[0]
    keep = np.delete(np.arange(0, m), delrow)
    A = A[keep, :]
    keep = np.delete(np.arange(0, m), delcol)
    A = A[:, keep]
    return A
 
# Max and min stiffness
Emin = 1e-9
Emax = 1.0
 
# Assemble FEM matrices and filter
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

# Degrees of Freedom matrix: edofMat[i,j,:] contains the
# DOFs associated to the square element (i,j)
    ndof = 2*(nelx+1)*(nely+1)
    dofs = np.arange(ndof).reshape((nelx+1, nely+1, 2))
    edofMat = np.zeros((nelx, nely, 8), dtype=int)
    edofMat[:, :, 0] = dofs[:-1, 1:, 0]
    edofMat[:, :, 1] = dofs[:-1, 1:, 1]
    edofMat[:, :, 2] = dofs[1:, 1:, 0]
    edofMat[:, :, 3] = dofs[1:, 1:, 1]
    edofMat[:, :, 4] = dofs[1:, :-1, 0]
    edofMat[:, :, 5] = dofs[1:, :-1, 1]
    edofMat[:, :, 6] = dofs[:-1, :-1, 0]
    edofMat[:, :, 7] = dofs[:-1, :-1, 1]
    edofMat = edofMat.reshape((nelx*nely, 8))

# FE: Build the index vectors for the for coo matrix format.
    KE = lk()
# Construct the index pointers for the coo format
    iK = np.kron(edofMat, np.ones((8, 1))).flatten()
    jK = np.kron(edofMat, np.ones((1, 8))).flatten()

# BC’s and support
    fixed = np.union1d(dofs[0, :, 0], dofs[-1,-1, 1])
    free = np.setdiff1d(dofs.reshape(ndof), fixed)

# Filter: Build (and assemble) the index+data vectors
# for the coo matrix format
    nfilter = int(nelx*nely*((2*(np.ceil(rmin)-1)+1)**2))
    iH = np.zeros(nfilter)
    jH = np.zeros(nfilter)
    sH = np.zeros(nfilter)
    cc = 0
    for i in range(nelx):
        for j in range(nely):
            row = i*nely+j
            kk1 = int(np.maximum(i-(np.ceil(rmin)-1), 0))
            kk2 = int(np.minimum(i+np.ceil(rmin), nelx))
            ll1 = int(np.maximum(j-(np.ceil(rmin)-1), 0))
            ll2 = int(np.minimum(j+np.ceil(rmin), nely))
            for k in range(kk1, kk2):
                for l in range(ll1, ll2):
                    col = k*nely+l
                    fac = rmin-np.sqrt(((i-k)*(i-k)+(j-l)*(j-l)))

                    iH[cc] = row
                    jH[cc] = col
                    sH[cc] = np.maximum(0.0, fac)
                    cc = cc+1
# Finalize assembly and convert to csc format
    H = sp.coo_matrix((sH, (iH, jH)),
        shape=(nelx*nely, nelx*nely)).tocsc()
    Hs = sp.diags(1/np.asarray(H.sum(1)).flatten())

# f is the unit load vector at the bottom right corner
    f = np.zeros((ndof, 1))
    f[1, 0] =-1

# Compute compliance and its sensitivity
@memoize()
def solve_state(x):
    dc = np.ones(nely*nelx)
    ce = np.ones(nely*nelx)

    # Setup and solve FE problem
    sK = ((KE.flatten()[np.newaxis]).T*(Emin+x**penal*(Emax-Emin)))\
    .flatten(order="F")
    K = sp.coo_matrix((sK, (iK, jK)), shape=(ndof, ndof)).tocsc()
    # Remove constrained dofs from matrix and convert to coo
    K = deleterowcol(K, fixed, fixed).tocoo()
    # Solve system
    K = cvxopt.spmatrix(K.data, K.row.astype(int), K.col.astype(int))
    B = cvxopt.matrix(f[free, :])
    cvxopt.cholmod.linsolve(K, B)
    u = np.zeros((ndof, f.shape[1]))
    u[free, :] = np.array(B)[:, :]

    # compliance
    ce = np.einsum("ija,jk,ika->ia", u[edofMat, :], KE, u[edofMat, :])
    obj = (Emin + (x**penal*(Emax-Emin))).dot(ce)

    # sensitivity of compliance with respect to x
    dc = np.multiply((-penal*x**(penal-1)*(Emax-Emin))[:, np.newaxis],
    ce).T

    # If not multi-objective
    if len(obj)==1:
        return (obj[0], dc[0,:])
    # else
    return (obj, dc)

# Density filter
@memoize()
def filter(x):
    return Hs @ (H @ x)


# Sensitivity of the filter
@memoize()
def diff_filter(x, v):
    return (H @ (Hs @ v.T)).T

# Definition of the optimization problem
@bound_constraints_optimizable(l=0, u=1)
@filtered_optimizable(filter, diff_filter)
class TO_problem(EuclideanOptimizable):
    def x0(self):
        return volfrac * np.ones(nely*nelx, dtype=float)

    def J(self, x):
        (obj, dc) = solve_state(x)
        return obj

    def dJ(self, x):
        (obj, dc) = solve_state(x)
        return dc

    def G(self, x):
        return [np.sum(x)/(nelx*nely)-volfrac]

    def dG(self, x):
        dv = np.ones(nelx*nely)/(nelx*nely)
        return dv

    def accept(self, params, results):
        # Plot the design at every iteration
        x = results["x"][-1]
        if not hasattr(self, "im"):
            plt.ion() # Ensure that redrawing is possible
            self.fig, self.ax = plt.subplots()
            self.im = self.ax.imshow(-x.reshape((nelx, nely)).T,\
            cmap="gray",
            interpolation="none",
            norm=colors.Normalize(vmin=-1, vmax=0))
            self.ax.axis("off")
            self.fig.show()
        else:
            self.im.set_array(-x.reshape((nelx, nely)).T)
            self.fig.canvas.draw()
        plt.pause(0.01)

# Optimization parameters
optimization_params = {"dt": 0.3,
                        "itnormalisation": 50,
                        "save_only_N_iterations": 1,
                        "save_only_Q_constraints": 5,

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