%% Multigrid test case with visulization
% Author: Ben Beaudry

clc; clear; close all
load A_b.mat;
pre = 2; % Number of presmoothing iterations
post = 2; % Number of postsmoothing iterations
% cycle = 1; % Type of multigrid cycle (1=V-cycle, 2=W-cycle, 3=F-cycle)
smooth = 1; % Smoother type (1=Jacobi, 2=Gauss-Seidel)
grids = 5; % Max grids in cycle
maxit = 1; % Max iterations of solver
tol = 1e-03; % Tolerance of solver

%% Solvers

% solve directly
t=cputime;
x_D = A\b;
fprintf('Direct Solve Time:        %g Seconds\n',cputime-t)

% V-Cycle
subplot(4,2,1);
hold on
title('V-Cycle')
t=cputime;
[x_V,res_V,it_V] = multigridVisual(A,b,pre,post,1,smooth,grids,maxit,tol);
fprintf('V-Cycle Solve Time:       %g Seconds\n',cputime-t)
grid on
hold off

% W-Cycle
subplot(4,2,[3,4]);
hold on
title('W-Cycle')
t=cputime;
[x_W,res_W,it_W] = multigridVisual(A,b,pre,post,2,smooth,grids,maxit,tol);
fprintf('W-Cycle Solve Time:       %g Seconds\n',cputime-t)
grid on
hold off

% F-Cycle
subplot(4,2,2);
hold on
title('F-Cycle')
t=cputime;
[x_F,res_F,it_F] = multigridVisual(A,b,pre,post,3,smooth,grids,maxit,tol);
fprintf('F-Cycle Solve Time:       %g Seconds\n',cputime-t)
grid on
hold off

% max it for iterative methods
it = 1000;

% Gauss-Siedel
t=cputime;
L = tril(A);
U = triu(A,1);
x_G = zeros(length(b),1);
res_G= zeros(it,1);
for g = 1:it
    x_G = L\(b-U*x_G);
    r = b - A * x_G;
    res_G(g) = norm(r);
    if res_G(g) < tol
        break;
    end
end
fprintf('Gauss-Seidel Solve Time:  %g Seconds\n',cputime-t)

% Jacobi
t=cputime;
d = diag(A);
dinv = 1./d;
LU = tril(A,-1)+triu(A,1);
U = triu(A,1);
x_J = zeros(length(b),1);
res_J= zeros(it,1);
for j = 1:it
    x_J = dinv.*(b-(LU*x_J));
    r = b - A * x_J;
    res_J(j) = norm(r);
    if res_J(j) < tol
        break;
    end
end
fprintf('Jacobi Solve Time:        %g Seconds\n',cputime-t)

%% Plotting

subplot(4,2,5:8);
hold on
plot(1:g,res_G(1:g),'o-','DisplayName','Guass-Seidel')
plot(1:j,res_J(1:j),'o-','DisplayName','Jacobi')
plot(1:it_V,res_V(1:it_V),'o-','DisplayName','V-Cycle Multigrid')
plot(1:it_F,res_F(1:it_F),'o-','DisplayName','F-Cycle Multigrid')
plot(1:it_W,res_W(1:it_W),'o-','DisplayName','W-Cycle Multigrid')
set(gca,'XLim',[0 100]);
grid on
legend('location','best')
ylabel('Relative Convergance')
xlabel('Iterations')
title('Convergance of Numerical System')
hold off