import numpy as np
from pyfreefem import FreeFemRunner
from pymedit import P0Function, P1Function
from nullspace_optimizer import EuclideanOptimizable, memoize, nlspace_solve,\
    bound_constraints_optimizable, filtered_optimizable
import matplotlib.pyplot as plt
import scipy.sparse.linalg as lg


def init(n=40, vfrac=0.1):
    global volfrac
    global Th
    global filter, diff_filter, solveA, MP1P0, solveMP0P0

    volfrac = vfrac

    mesh_script = """   
        IMPORT "io.edp"
        mesh Th = square($n, $n, flags=1);  
        exportMesh(Th);"""
    Th = FreeFemRunner(mesh_script).execute({'n': n})['Th']

    # Changing boundary label on left hand side
    pts1 = Th.vertices[Th.edges[:, 0]-1]
    pts2 = Th.vertices[Th.edges[:, 1]-1]
    eps = 1/(10*n)
    indices = np.logical_and(
        abs(pts1[:, 1]-0.5) <= 0.1+eps, abs(pts2[:, 1]-0.5) <= 0.1+eps)
    indices = np.logical_and(indices, Th.edges[:, -1] == 4)
    Th.edges[indices, -1] = 5
    Th.edges[Th.edges[:, -1] != 5, -1] = 1
    Th._AbstractMesh__updateBoundaries()
    # Th.plot(bcEdgeLabels=True)

    filter_script = """ 
    IMPORT "io.edp" 
    mesh Th = importMesh("Th");
            
    fespace Fh0(Th,P0); 
    fespace Fh1(Th,P1);     

    varf mass(u,v)=int2d(Th)(u*v); 

    matrix MP1P0 = mass(Fh0,Fh1);   
    matrix MP0P0 = mass(Fh0,Fh0);   
        
    macro grad(u) [dx(u),dy(u)]//   
            
    real gamma = Th.hmin; 
    varf helmholtz(u,v) = int2d(Th)(gamma^2*grad(u)'*grad(v)+u*v);    
    matrix A = helmholtz(Fh1,Fh1);  
        
    exportMatrix(MP1P0);    
    exportMatrix(MP0P0);    
    exportMatrix(A); 
    """

    runner = FreeFemRunner(filter_script)
    runner.import_variables(Th=Th)
    exports = runner.execute()
    solveA = lg.factorized(exports['A'])
    solveMP0P0 = lg.factorized(exports['MP0P0'])
    MP1P0 = exports['MP1P0']


init(5)


@memoize()
def filter(rho):
    return solveMP0P0(MP1P0.T @ solveA(MP1P0 @ rho))


@memoize()
def diff_filter(rho, v):
    return (MP1P0.T @ solveA(MP1P0 @ solveMP0P0(v.T))).T


@memoize()
def solve_state(rho):
    solve_script = """  
    IMPORT "io.edp" 
        
    mesh Th = importMesh("Th");
        
    fespace Fh0(Th,P0); 
    fespace Fh1(Th,P1);     
        
    Fh0 rho;    
    rho[] = importArray("rho");
        
    real p = $p;
    real kappaf = $kappaf;  
    real kappas = $kappas;
    Fh0 kappa = rho^p*(kappaf-kappas)+kappas;
    Fh0 dkappa = p*rho^(p-1)*(kappaf-kappas);
    Fh1 T, S;   
    func Q = $Q;
        
    macro grad(u) [dx(u),dy(u)] //
        
    solve heat(T,S)=    
        int2d(Th)(kappa*grad(T)'*grad(S))   
        -int2d(Th)(Q*S)
        +on(5,T=0);
            
    real J = int2d(Th)(T); 
    real vol0 = int2d(Th)(1.);
    real volfrac = $volfrac;
    real H = int2d(Th)(rho/(vol0*volfrac))-1;  
        
    Fh0 drho, dummy; 
    varf DJ(dummy, drho) = -int2d(Th)(dkappa*grad(T)'*grad(T)*drho*1.0/Q);
    real[int] dJ = DJ(0,Fh0);

    varf DH(dummy, drho) = int2d(Th)(drho/(vol0*volfrac));
    real[int] dH = DH(0,Fh0);

    exportVar(J); 
    exportVar(H);
    exportArray(T[]);
    exportArray(dJ);
    exportArray(dH);
    """
    runner = FreeFemRunner(solve_script, debug=0)
    runner.import_variables(rho=rho, Th=Th)
    exports = runner.execute({'p': 3, 'Q': 1e4,
                              'kappaf': 401,
                              'kappas': 1,
                              'volfrac': volfrac})
    return exports


@bound_constraints_optimizable(l=0, u=1)
@filtered_optimizable(filter, diff_filter)
class Heat_TO(EuclideanOptimizable):
    def __init__(self, plot=None):
        self.plot = plot

        if self.plot:
            plt.ion()
            self.fig, self.ax = plt.subplots(1, 2)
            self.fig.set_size_inches(12, 6)
            self.ax = self.ax.flatten()
            self.cbars = [None]*2

    def x0(self):
        return volfrac * np.ones(Th.nt, dtype=float)

    def J(self, rho):
        return solve_state(rho)['J']

    def H(self, rho):
        return [solve_state(rho)['H']]

    def dJ(self, rho):
        return solve_state(rho)['dJ']

    def dH(self, rho):
        return solve_state(rho)['dH']

    def accept(self, params, results):
        if self.plot:
            #self.ax[0].clear()

            rho = results['x'][-1]
            for i in range(2):
                if self.cbars[i]:
                    self.cbars[i].remove()
                self.ax[i].clear()
            _, _, self.cbars[0] = P0Function(Th, rho).plot(fig=self.fig, ax=self.ax[0], cmap="gray_r",
                                                           vmin=0, vmax=1)
            _, _, self.cbars[1] = P1Function(Th, solve_state(rho)['T[]']).plot(fig=self.fig, cmap="turbo",
                                                                               ax=self.ax[1],
                                                                               type_plot='tricontourf')
            self.ax[0].title.set_text('rho')
            self.ax[1].title.set_text('T')
            self.fig.suptitle('Iteration '+str(results['it'][-1]))
            # savefig(fig_dir+"/"+it+".png", fig=self.fig, level=1, sizefactor=1)
            plt.pause(0.1)


if __name__ == "__main__":
    init(100)
    plt.ion()

    case = Heat_TO(plot=True)
    params = dict(dt=0.1, itnormalisation=50, maxit=150,
                  save_only_N_iterations=1,
                  save_only_Q_constraints=5,
                  qp_solver='qpalm',
                  show_progress_qp=False)
    nlspace_solve(case, params)
