%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  lab  Iterations methods  %%%
%%%   for elliptic equation   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

% task params
center = [0, 0];
edge_len = 2;

% coefficients of equation
a_coeff = @(x, y) ones(size(y))'*x.^2;
b_coeff = @(x, y) y.^2'*ones(size(x));
q_coeff = @(x, y) 60*ones(size(y))'*ones(size(x)) + y.^2'*x.^2;
% right side of equation
f = @(x, y) y.^7'*x.^7;

% exact solution
u_ex = @(x, y) y.^5'*x.^5;
% border conditions
border = u_ex;

% task discretization:
level = 8;
% number of inner nodes by one direction
N = 2^level-1;

% init solver
Solver = ellipticProblemSolver(center, edge_len, N,...
                               a_coeff, b_coeff, q_coeff,...
                               f, border, u_ex);

% set solver params
eps = 1e-4;
k_max = 20;
Solver.setParams(eps, k_max);

% solve by miltigrid method
nu = 3;
level = 3;
smooth_ms = {'Jacobi', 'Seidel', 'SOR'};
v = Solver.solveProblemByMG(nu, nu, level, char(smooth_ms{1}));
       
% plot results
figure
subplot(2, 1, 1)
indx = 1:length(Solver.error);
semilogy(indx, Solver.error, 'o-')
grid on
title('Error')

subplot(2, 1, 2)
plot(indx, Solver.rho, '*-')
grid on
title('Spectral radius')

