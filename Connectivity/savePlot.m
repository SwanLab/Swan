% yyaxis left
plot(cp1.cost); hold on;  
plot(cp08.cost);  
plot(cp06.cost);
% yyaxis right
% yyaxis right
plot(-p1.constraint(:,2)+1); hold on;  
plot(-p08.constraint(:,2)+0.8);  
plot(-p06.constraint(:,2)+0.6);
ylabel('First Eigenvalue', 'FontSize', 14)
xlabel('Iteration', 'FontSize', 14)
grid on
legend({'$\lambda_1^{\min} = 1.0$ ','$\lambda_1^{\min} = 0.8$ ','$\lambda_1^{\min} = 0.6$ '},'Interpreter','latex','Location','northeast', 'FontSize', 14)

figure; 
plot(-np1.constraint(:,2)+1); hold on;  
plot(-np08.constraint(:,2)+0.8);  
plot(-np06.constraint(:,2)+0.6);
ylabel('First Eigenvalue', 'FontSize', 14)
xlabel('Iteration', 'FontSize', 14)
grid on
legend({'$\lambda_1^{\min} = 1.0$ ','$\lambda_1^{\min} = 0.8$ ','$\lambda_1^{\min} = 0.6$ '},'Interpreter','latex','Location','northeast', 'FontSize', 14)



% figure;
% iterations = 1:1:600
% yyaxis left
% plot(iterations, -cost(1:600), '-')
% ylabel('First Eigenvalue', 'FontSize', 14)
% yyaxis right
% 
% plot(iterations,  constraint(1:600), '-.');hold on
% % semilogy(iterations,  max, '--')
% xlabel('Iteration', 'FontSize', 14)
% grid
% legend({'$\lambda_1$','$\|T\|_{2}$','$\|T\|_{\infty}$'},'Interpreter','latex','Location','southeast', 'FontSize', 14)
% ax = gca;          % get current axes
% ax.FontSize = 12;

% load('cost.mat')
% load('constraint.mat')
% % 
% per = load('PDEincreasingConstraintPerimeter.mat');
% noConst= load('PDEnoConst.mat');
% cons = load('PDEincreasingConstraint.mat');
% 
% cost1 = load('cost_cantilever_case_1.mat')
% constraint1 = load('constraint_cantilever_case_1.mat')
% % cost2 = load('cost_cantilever_case_2.mat')
% % constraint2 = load('constraint_cantilever_case_2.mat')
% cost3 = load('cost_cantilever_case_3.mat')
% constraint3 = load('constraint_cantilever_case_3.mat')
% cost4 = load('cost_cantilever_case_4.mat')
% constraint4 = load('constraint_cantilever_case_4.mat')
% cost5 = load('cost_cantilever_case_5.mat')
% constraint5 = load('constraint_cantilever_case_5.mat')
% 
% figure;
% N = 1000;
% iterations = 1:1:N;
% % lambda1min = min((0.5 * kron(0:ceil(N/20)-1, ones(1,20))), 6.0); hold on
% % lambda1min = min(0.4+(0.2 * kron(0:ceil(N/20)-1, ones(1,20))), 4.0)
% lambda1min = min(0.4+(1.0 * kron(0:ceil(N/20)-1, ones(1,20))), 3.4);

% figure
% semilogy(iterations, -noCost.cost, 'k-'); hold on
% semilogy(iterations, 2.0 - constraintA.constraint(1:1000,2), '-'); hold on
% semilogy(iterations, 1.0 - constraintB.constraint(1:1000,2), '-');
% semilogy(iterations, 2.0 - constraintC.constraint(1:1000,2), '-');
% semilogy(iterations, 2.0 - constraintD.constraint(1:1000,2), '-');
% semilogy(iterations, lambda1min(:) - constraintE.constraint(1:1000,2), '-');

% semilogy(iterations, lambda1min(:) - constraint(:,2), '-')
% semilogy(iterations, -noConst.cost, 'k-'); hold on
% semilogy(iterations, 0.4 - fixed04.constraint(1:500,2), '-')
% semilogy(iterations, 4 - fixed4.constraint(1:500,2), '-')
% semilogy(iterations, lambda1min(:) - increasing.constraint(1:500,2), '-')
% semilogy(iterations, 0.4*ones(size(iterations)),'--')
% semilogy(iterations, 4.0*ones(size(iterations)),'--')
% lambda1min = min(0.2+(0.2 * kron(0:ceil(N/20)-1, ones(1,20))), 4.0)
% semilogy(iterations, lambda1min(:) - perIncreasing.constraint(1:500,2), '-')

% 
% ylabel('First Eigenvalue')
% xlabel('Iteration')
% % legend({'No connectivity constraints','$\lambda_1^{\min} = 0.4$ fixed','$\lambda_1^{\min} = 4.0$ fixed','$\lambda_1^{\min}$ increasing from 0.2 to 4.0','$\lambda_1^{\min}$ increasing from 0.2 to 4.0 + perimeter penalization'},'Location','southoutside','Interpreter','latex')
% legend({'No connectivity constraints','Mass + $\lambda_1^{\min}=2.0$ fixed','PDE + $\lambda_1^{\min}=1.0$ fixed','PDE + $\lambda_1^{\min}=2.0$ fixed + perimeter penalization','PDE + HP + $\lambda_1^{\min}=2.0$ fixed + perimeter penalization','PDE + HP + $\lambda_1^{\min}$ increasing from 0.2 to 4.0 + perimeter penalization'},'Location','southoutside','Interpreter','latex')
% grid

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
