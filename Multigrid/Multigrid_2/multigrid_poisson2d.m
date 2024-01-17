

function u  = multigrid_poisson2d(u, f, L, Ncycles, Npre,N)
% Inputs"
% N - number of grid points in each direction (minus boundary points)
% nu1 - number of pre-smoothing iterations
% nu1 - number of post-smoothing iterations
% levels - number of levels in the V-cycle
% v1 - initial guess for the solution
% v2 - right hand side vector of the linear system

% Outputs:
% u - approximate solution to the linear system
% res - vector of residual norms at each iteration

%