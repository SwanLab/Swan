Colors = get(groot, 'defaultAxesColorOrder');
% set(groot,'DefaultAxesFontSize',14);
% set(groot,'DefaultTextFontSize',15);
% set(groot,'DefaultLineLineWidth',1.2);
% set(groot,'DefaultAxesLineWidth',1.2);

set(groot,'DefaultAxesFontSize',14);     % axis numbers
set(groot,'DefaultAxesLabelFontSizeMultiplier',1.2);
set(groot,'DefaultAxesTitleFontSizeMultiplier',1.3);
set(groot,'DefaultTextFontSize',14);     % text, legends, annotations
set(groot,'DefaultLegendFontSize',14);
set(groot,'DefaultLineLineWidth',1.2);   % thicker lines


% %%%%% ONLY COMPLIANCE, VOLUME, EIG %%%%%%%%%%%%%
fig = figure('Units','centimeters','Position',[1 2 24 6]);
t = tiledlayout(1,2);
t.TileSpacing = 'loose';
t.Padding = 'compact';
% 
% % =====================================================
% % TOP LEFT: Compliance
% =====================================================
nexttile
plot(-1.*cost(1:3000),  'Color', Colors(1,:), 'LineStyle','-'); hold on
% plot((constraint(1:2000,1)+1)*0.4, 'Color', Colors(2,:), 'LineStyle',':');
grid on
xlabel('Iteration'); 
ylabel('First Eigenvalue'); 
ylim([-1 80])
% legend({'First Eigenvalue','Volume'}, 'Location','best','Box','off');

% =====================================================
% BOTTOM CENTERED: First Eigenvalue
% =====================================================
nexttile
plot((constraint(1:3000,1)+1)*0.4,  'Color', Colors(2,:), 'LineStyle','-'); hold on;ylim([0,1])
% plot(constraint2.constraint(:,1), 'Color', Colors(2,:), 'LineStyle','-');
grid on
xlabel('Iteration'); ylabel('Volume');
% legend({'Mass','PDE'}, 'Location','best','Box','off');

% =====================================================
% % BOTTOM CENTERED: First Eigenvalue
% % =====================================================
% ax3 = axes('Parent',fig,'Position',[xCenter yBottom w h]);
% plot(-constraint(:,2)+1,  'Color', Colors(5,:), 'LineStyle','-'); hold on
% % plot(-constraint2.constraint(:,2)+1, 'Color', Colors(5,:), 'LineStyle','-');
% grid on
% xlabel('Iteration'); ylabel('First Eigenvalue');
% % legend({'Mass','PDE'}, 'Location','best','Box','off');

% =====================================================
% EXPORT PDF WITHOUT MARGINS
% =====================================================
figPos = get(fig,'Position');
set(fig,'PaperUnits','centimeters');
set(fig,'PaperSize',figPos(3:4));
set(fig,'PaperPosition',[0 0 figPos(3:4)]);

print(fig, 'eigMax.pdf', '-dpdf', '-painters');

% % 
% % %%%%% ONLY COMPLIANCE, VOLUME, EIG %%%%%%%%%%%%%
% fig = figure('Units','centimeters','Position',[2 2 24 12]);
% % ---- layout parameters (normalized) ----
% w  = 0.38;           % width of each plot  (same for all)
% h  = 0.32;           % height of each plot (same for all)
% 
% xLeft  = 0.07;       % x for plot 1 (top-left)
% xRight = 1 - xLeft - w;   % x for plot 2 (top-right)
% 
% % center position for bottom plot:
% xCenter = 0.5 - w/2;   
% 
% yTop    = 0.57;      % y for top row
% yBottom = 0.12;      % y for bottom row
% 
% % =====================================================
% % TOP LEFT: Compliance
% % =====================================================
% ax1 = axes('Parent',fig,'Position',[xLeft yTop w h]);
% plot(costMass(1:1000),'k:'); hold on
% plot(costPDE(1:1000), 'Color', Colors(1,:), 'LineStyle','-');
% grid on
% xlabel('Iteration'); ylabel('Compliance');
% legend({'Mass','PDE'}, 'Location','best','Box','off');
% 
% % =====================================================
% % BOTTOM CENTERED: First Eigenvalue
% % =====================================================
% ax2 = axes('Parent',fig,'Position',[xRight yTop w h]);
% plot((constraintMass(1:1000,1)+1)*0.4, 'k:' ); hold on
% plot((constraintPDE(1:1000,1)+1)*0.4,'Color', Colors(2,:), 'LineStyle','-' ); ylim([0,1])
% grid on
% xlabel('Iteration'); ylabel('Volume');
% legend({'Mass','PDE'}, 'Location','best','Box','off');
% 
% % =====================================================
% % BOTTOM CENTERED: First Eigenvalue
% % =====================================================
% ax3 = axes('Parent',fig,'Position',[xCenter yBottom w h]);
% plot(-constraintMass(1:100,2)+1, 'k:'); hold on
% plot(-constraintPDE(1:1000,2)+1, 'Color', Colors(5,:),  'LineStyle','-');
% grid on
% xlabel('Iteration'); ylabel('First Eigenvalue');
% legend({'Mass','PDE'}, 'Location','best','Box','off');
% 
% % =====================================================
% % EXPORT PDF WITHOUT MARGINS
% % =====================================================
% figPos = get(fig,'Position');
% set(fig,'PaperUnits','centimeters');
% set(fig,'PaperSize',figPos(3:4));
% set(fig,'PaperPosition',[0 0 figPos(3:4)]);
% 
% print(fig, 'massPDE.pdf', '-dpdf', '-painters');


% %%%%% COMPLIANCE, PERIMETER, VOLUME, EIG %%%%%%%%%%%%%
% % 
% fig = figure('Units','centimeters','Position',[2 2 24 12]);
% t = tiledlayout(2,2);
% t.TileSpacing = 'loose';
% t.Padding = 'compact';
% 
% % ----- TOP LEFT: Compliance -----
% nexttile
% plot(cost,  'k:');  hold on
% plot(costPDE, 'Color', Colors(1,:),'LineStyle','-');
% grid on
% xlabel('Iteration'); ylabel('Compliance');
% legend({'Mass','PDE'}, 'Location','best', 'Box','off');
% 
% % ----- TOP RIGHT: Perimeter -----
% nexttile
% plot(per(1:1000),  'k:'); hold on
% plot(perPDE(1:1000), 'Color', Colors(3,:), 'LineStyle','-');
% grid on
% xlabel('Iteration'); ylabel('Perimeter');
% legend({'Mass','PDE'}, 'Location','best', 'Box','off');
% 
% % ----- BOTTOM LEFT: Volume -----
% nexttile
% plot((constraint(:,1)+1)*0.4,   'k:'); hold on
% plot((constPDE(:,1)+1)*0.4,  'Color', Colors(2,:),  'LineStyle','-'); ylim([0.0,1.0])
% grid on
% xlabel('Iteration'); ylabel('Volume');
% legend({'Mass','PDE'}, 'Location','best', 'Box','off');
% 
% % ----- BOTTOM RIGHT: First Eigenvalue -----
% nexttile  % <---- THIS CENTERS PLOT #3
% plot(-constraint(:,2)+1,   'k:'); hold on
% plot(-constPDE(:,2)+1,  'Color', Colors(5,:), 'LineWidth',1.6, 'LineStyle','-');
% grid on
% xlabel('Iteration'); ylabel('First Eigenvalue');
% legend({'Mass','PDE'}, 'Location','best', 'Box','off');
% 
% 
% % ----- EXPORT PDF WITHOUT MARGINS -----
% fig = gcf;
% figPos = get(fig,'Position');          % [x y width height]
% set(fig,'PaperUnits','centimeters');
% set(fig,'PaperSize',figPos(3:4));
% set(fig,'PaperPosition',[0 0 figPos(3:4)]);
% 
% print(fig, 'Perimeter.pdf', '-dpdf', '-painters');

% print(figure(1),'figure2x2','-dpdf');   % best for papers
% exportgraphics(gcf, 'myFigure.pdf', 'ContentType','vector', 'BackgroundColor','none');

% set(groot, 'DefaultAxesFontSize', 16);        % axis numbers
% set(groot, 'DefaultAxesLabelFontSizeMultiplier', 1.2);
% set(groot, 'DefaultAxesTitleFontSizeMultiplier', 1.2);
% 
% set(groot, 'DefaultLineLineWidth', 1.8);      % all lines
% set(groot, 'DefaultAxesLineWidth', 1.2);      % axis box
% 
% set(groot, 'DefaultTextFontSize', 16);        % text, legend, colorbar
% set(groot, 'DefaultLegendFontSize', 16);
% 
% set(groot, 'DefaultAxesFontName', 'Times');   % optional (for papers)
% set(groot, 'DefaultTextFontName', 'Times');
% 
% figure(1);
% tiledlayout(2,2)
% nexttileti
% plot(cost,'Color', Colors(1,:)); ylabel('Compliance'); xlabel('Iteration'); grid on
% nexttile
% plot(cost,'Color', Colors(2,:)); ylabel('Perimeter'); xlabel('Iteration'); grid on
% 
% nexttile
% plot(constraint(:,1),'Color', Colors(3,:)); ylabel('Volume'); xlabel('Iteration'); grid on
% nexttile
% plot(-constraint(:,2)+1.0,'Color', Colors(4,:)); ylabel('First Eigenvalue'); xlabel('Iteration'); grid on
% saveas(figure(1),'PDEperimeter','png')
% 

% % ---- Figure 2 ----
% f2 = figure(2);
% % yyaxis left
% plot(constraint(:,1), 'Color', Colors(5,:)); hold on;
% ylabel('Volume', 'FontSize', 14)
% 
% % yyaxis right
% plot(-constraint(:,2)+1, 'Color', Colors(4,:)); hold on;
% ylabel('First Eigenvalue', 'FontSize', 14)
% xlabel('Iteration', 'FontSize', 14)
% grid on
% 
% % ax = gca;
% % ax.YAxis(1).Color = Colors(5,:);  % left axis same as left line
% % ax.YAxis(2).Color = Colors(4,:);  % right axis same as right line
% 
% saveas(f2,'constraintPDEperimeter','png')
% 

% % yyaxis left
% plot(cp1.cost); hold on;  
% plot(cp08.cost);  
% plot(cp06.cost);
% % yyaxis right
% % yyaxis right
% plot(-p1.constraint(:,2)+1); hold on;  
% plot(-p08.constraint(:,2)+0.8);  
% plot(-p06.constraint(:,2)+0.6);
% ylabel('First Eigenvalue', 'FontSize', 14)
% xlabel('Iteration', 'FontSize', 14)
% grid on
% legend({'$\lambda_1^{\min} = 1.0$ ','$\lambda_1^{\min} = 0.8$ ','$\lambda_1^{\min} = 0.6$ '},'Interpreter','latex','Location','northeast', 'FontSize', 14)
% 
% figure; 
% plot(-np1.constraint(:,2)+1); hold on;  
% plot(-np08.constraint(:,2)+0.8);  
% plot(-np06.constraint(:,2)+0.6);
% ylabel('First Eigenvalue', 'FontSize', 14)
% xlabel('Iteration', 'FontSize', 14)
% grid on
% legend({'$\lambda_1^{\min} = 1.0$ ','$\lambda_1^{\min} = 0.8$ ','$\lambda_1^{\min} = 0.6$ '},'Interpreter','latex','Location','northeast', 'FontSize', 14)



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
