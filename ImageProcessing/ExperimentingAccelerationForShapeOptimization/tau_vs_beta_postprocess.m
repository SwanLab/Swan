%% POSTPROCESS
clear
% close all
load("Tau_vs_Beta_bridge1.mat");
iters_tau_beta_constant = iters;
iters_tau_beta_constant(iters_tau_beta_constant<20) = 60;
load("adaptative_beta_bridge1.mat");
iters_adaptative = iters;
figure()
hold on
s = surf(betaGrid,tauGrid,iters_tau_beta_constant);
xlabel('Momentum term ($\beta$)','Interpreter','latex')
ylabel('Line search ($\tau$)','Interpreter','latex')
s.EdgeColor = "none";
s.FaceColor = "interp";
c = colorbar;
colormap(flipud(jet));
c.TickLabelInterpreter = 'latex';
clim([min(min(iters_tau_beta_constant)),60]);
view(0,90)
bestOnes = zeros(size(iters_tau_beta_constant,1),3);
for j = 1:size(iters_tau_beta_constant,1)
    [val,pos] = min(iters_tau_beta_constant(j,:));
    bestOnes(j,:) = [betaGrid(1,pos),tauGrid(j,1),val+1e3];
end
h = plot3(bestOnes(:,1),bestOnes(:,2),bestOnes(:,3),'o','Color','black');
set(h, 'MarkerFaceColor', 'k'); 
set(gca,"TickLabelInterpreter",'latex','FontSize',14)
xlim([0 1])
ylim([50 200])
box on
hold off
figure()
hold on
plot(bestOnes(:,2),bestOnes(:,3)-1e3,'-o','Color','black')
plot(tau,iters_adaptative,'--*','Color','black')
xlabel('Line search ($\tau$)','Interpreter','latex')
ylabel('Iterations to converge','Interpreter','latex')
legend('Optimum constant $\beta$','Adaptative $\beta$','Interpreter','latex')
set(gca,"TickLabelInterpreter",'latex','FontSize',14)
ylim([0,100])
xlim([50 200])
grid minor
box on
hold off
