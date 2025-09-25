% load('cost.mat')
% load('constraint.mat')
% % 
% per = load('PDEincreasingConstraintPerimeter.mat');
% noConst= load('PDEnoConst.mat');
% cons = load('PDEincreasingConstraint.mat');

figure;
N = 500;
iterations = 1:1:N;
% lambda1min = min((0.5 * kron(0:ceil(N/20)-1, ones(1,20))), 6.0); hold on
lambda1min = min(0.4+(0.2 * kron(0:ceil(N/20)-1, ones(1,20))), 4.0)

% semilogy(iterations, lambda1min(:) - constraint(:,2), '-')
semilogy(iterations, -noConst.cost, 'k-'); hold on
semilogy(iterations, 0.4 - fixed04.constraint(1:500,2), '-')
semilogy(iterations, 4 - fixed4.constraint(1:500,2), '-')
semilogy(iterations, lambda1min(:) - increasing.constraint(1:500,2), '-')
% semilogy(iterations, 0.4*ones(size(iterations)),'--')
% semilogy(iterations, 4.0*ones(size(iterations)),'--')

lambda1min = min(0.2+(0.2 * kron(0:ceil(N/20)-1, ones(1,20))), 4.0)
semilogy(iterations, lambda1min(:) - perIncreasing.constraint(1:500,2), '-')


ylabel('First Eigenvalue')
xlabel('Iteration')
legend({'No connectivity constraints','$\lambda_1^{\min} = 0.4$ fixed','$\lambda_1^{\min} = 4.0$ fixed','$\lambda_1^{\min}$ increasing from 0.2 to 4.0','$\lambda_1^{\min}$ increasing from 0.2 to 4.0 + perimeter penalization'},'Location','southoutside','Interpreter','latex')
grid

% save('PDEincreasingCost.mat','cost')
% save('PDEincreasingConstraint.mat','constraint')

% figure;
% iterations = 1:1:500
% yyaxis left
% semilogy(iterations, -cost, '-.')
% ylabel('First Eigenvalue')
% yyaxis right
% semilogy(iterations, 0.4*(constraint+1), '-.')
% ylabel('Constraint')
% xlabel('Iteration')
% grid
% xlabel('Iterations')

% load('PDE.mat'); PDE = lambda1
% load('LUMP.mat'); LUMP = lambda1
% load('PDEFP.mat'); PDEFP = lambda1
% 
% figure;
% 
% h = 1/100
% semilogy(0.0:h:1.0, LUMP, 'r-'); hold on
% semilogy(0.0:h:1.0, PDE, 'b-'); hold on
% semilogy(0.0:h:1.0, PDEFP, 'g-'); hold on
% ylabel('First Eigenvalue')
% xlabel('Radius')
% grid
