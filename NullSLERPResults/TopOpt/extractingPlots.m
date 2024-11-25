
clear;
close all;

% 3 gJs case
fig1 = openfig('NullSLERPResults/TopOpt/2DCantileverBeam/FinalResults_NoOscillations/Monitoring_trust0d02_gJ0d2.fig');
fig2 = openfig('NullSLERPResults/TopOpt/2DCantileverBeam/FinalResults_NoOscillations/Monitoring_trust0d02_gJ1.fig');
fig3 = openfig('NullSLERPResults/TopOpt/2DCantileverBeam/FinalResults_NoOscillations/Monitoring_trust0d02_gJ2d5V0d4.fig');

iter  = fig1.Children(end).Children.XData;
cost1 = fig1.Children(end).Children.YData;
cost2 = fig2.Children(end).Children.YData;
cost3 = fig3.Children(end).Children.YData;

constr1 = fig1.Children(end-2).Children.YData;
constr2 = fig2.Children(end-2).Children.YData;
constr3 = fig3.Children(end-2).Children.YData;

lambda1 = fig1.Children(end-4).Children.YData;
lambda2 = fig2.Children(end-4).Children.YData;
lambda3 = fig3.Children(end-4).Children.YData;

eta1 = fig1.Children(end-7).Children.YData;
eta2 = fig2.Children(end-7).Children.YData;
eta3 = fig3.Children(end-7).Children.YData;

close all;

colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]};

f1 = figure;
plot(iter,cost1,iter,cost2,iter,cost3,'LineWidth',1.2);
grid on
grid minor
legend('$\eta^*=0.20$','$\eta^*=1.00$','$\eta^*=2.50$','interpreter','latex','FontSize', 12.5)
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f1.Position = [680   558   560/1.7   420];

f2 = figure;
plot(iter,constr1,iter,constr2,iter,constr3,'LineWidth',1.2)
grid on
grid minor
legend('$\eta^*=0.20$','$\eta^*=1.00$','$\eta^*=2.50$','interpreter','latex','FontSize', 12.5)
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f2.Position = [680   558   560/1.7   420];

f3 = figure;
plot(iter,lambda1,iter,lambda2,iter,lambda3,'LineWidth',1.2)
grid on
grid minor
legend('$\eta^*=0.20$','$\eta^*=1.00$','$\eta^*=2.50$','interpreter','latex','FontSize', 12.5)
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f3.Position = [680   558   560/1.7   420];

f4 = figure;
h = plot(iter(1:162),eta1(1:162),'--',iter(163:end),eta1(163:end),iter(1:74),eta2(1:74),'--',iter(75:end),eta2(75:end),iter(1:50),eta3(1:50),'--',iter(51:end),eta3(51:end),'LineWidth',1.2);
grid on
grid minor
legend('','$\eta^*=0.20$','','$\eta^*=1.00$','','$\eta^*=2.50$','interpreter','latex','FontSize', 12.5)
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f4.Position = [680   558   560/1.7   420];
[h(1).Color, h(3).Color, h(5).Color] = colors{:};
[h(2).Color, h(4).Color, h(6).Color] = colors{:};