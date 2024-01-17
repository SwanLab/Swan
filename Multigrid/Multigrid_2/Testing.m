close all 
clear all
%Define problem parameterts

L = 129; % number of grid points in each direction
[X, Y] = meshgrid(linspace(0, 1, L), linspace(0, 1, L));

f = rhs_func(X,Y);
u_exact = exact_soultion(X,Y);

%Set up the initial guess


%disp(size(uo))

%Set up multigrid parameters
Ncycles = 11;
Npre = 23;
Npost = 23;

%disp(size(f))
%disp(u_exact)
%disp(size(u_exact))
uo = zeros(L,L);

% Solve using multigrid V -cycle with Jacobi smoother and bilinear
% interpolation

%u = multigrid_v_cycle_jacobi(uo, f, L, Ncycles, Npre, Npost);
U{1} = multigrid_v_cycle_jacobi(uo, f, L, Ncycles, Npre, Npost);
%disp(U{1})
% Plot numerical solution and exact solution
 
figure;
% subplot(1, 2, 1);
% surf(X, Y, u);
% title('Numerical solution');
% xlabel('x');
% ylabel('y');
% zlabel('u');
% subplot(1, 2, 2);
surf(X, Y, u_exact);
title('Exact solution');
xlabel('x');
ylabel('y');
zlabel('u');
% 
