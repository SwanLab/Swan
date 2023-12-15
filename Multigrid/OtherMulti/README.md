# Multigrid method for elliptic equations

This study project includes the entire algorithm of the [Multigrid method](https://en.wikipedia.org/wiki/Multigrid_method), which was applied to simple math problem as an example.
The problem contains an elliptic differential equation, Dirichlet boundary condition on a square domain; exact solution is known.
The problem is approximated by the finite volume method.
One of three smoothing methods can be chosen: Jacobi, Seidel, SOR.
On the coarsest grid the problem is solved by block tridiagonal method.
The spectral radius of iterative matrix is calculated during the solving.
