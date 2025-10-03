from nullspace_optimizer.optimizers.MMA import mma_lib, mma_solve
from nullspace_optimizer.examples.topopt_examples import ex00_compliance as ex00
from nullspace_optimizer.examples.topopt_examples import ex03_compliance_MMA as ex03

ex00.init(nelx=150, nely=50, volfrac=0.01, rmin=1.9)

if __name__ == "__main__":
    case = ex03.TO_problem(plot=True)
    mma_lib.FIX = True  # Fix initialization of the asymptotes as per Niels Aage
    mma_solve(case, 0, 1, {'move': 0.2,
                           'save_only_N_iterations': 5,
                           'maxit': 150,
                           'tight_move': True})
