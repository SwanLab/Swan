%% 3 gJs case

clear;
close all;

fig1 = openfig('NullSLERPResults/TopOpt/MBBBeam/FinalResults_NoOscillations/Monitoring_trust0d02_gJ0.5V0d4.fig');
fig2 = openfig('NullSLERPResults/TopOpt/MBBBeam/FinalResults_NoOscillations/Monitoring_trust0d02_gJ1V0d4.fig');
fig3 = openfig('NullSLERPResults/TopOpt/MBBBeam/FinalResults_NoOscillations/Monitoring_trust0d02_gJ2V0d4.fig');

iter1 = fig1.Children(end).Children.XData;
iter2 = fig2.Children(end).Children.XData;
iter3 = fig3.Children(end).Children.XData;

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
plot(iter1,cost1,iter2,cost2,iter3,cost3,'LineWidth',1.2);
grid on
grid minor
legend('$\eta^*=0.50$','$\eta^*=1.00$','$\eta^*=2.00$','interpreter','latex','FontSize', 12.5)
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f1.Position = [680   558   560/1.7   420];

f2 = figure;
plot(iter1,constr1,iter2,constr2,iter3,constr3,'LineWidth',1.2)
grid on
grid minor
legend('$\eta^*=0.50$','$\eta^*=1.00$','$\eta^*=2.00$','interpreter','latex','FontSize', 12.5)
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f2.Position = [680   558   560/1.7   420];

f3 = figure;
plot(iter1,lambda1,iter2,lambda2,iter3,lambda3,'LineWidth',1.2)
grid on
grid minor
legend('$\eta^*=0.50$','$\eta^*=1.00$','$\eta^*=2.00$','interpreter','latex','FontSize', 12.5)
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f3.Position = [680   558   560/1.7   420];

f4 = figure;
h = plot(iter1(1:53),eta1(1:53),'--',iter1(54:end),eta1(54:end),iter2(1:48),eta2(1:48),'--',iter2(49:end),eta2(49:end),iter3(1:57),eta3(1:57),'--',iter3(58:end),eta3(58:end),'LineWidth',1.2);
grid on
grid minor
legend('','$\eta^*=0.50$','','$\eta^*=1.00$','','$\eta^*=2.00$','interpreter','latex','FontSize', 12.5)
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f4.Position = [680   558   560/1.7   420];
[h(1).Color, h(3).Color, h(5).Color] = colors{:};
[h(2).Color, h(4).Color, h(6).Color] = colors{:};


%% 2 gJs case

clear;
close all;

fig1 = openfig('NullSLERPResults/TopOpt/Gripper/FinalResults_NoOscillations/Monitoring_trust0d02_gJ0.05.fig');
fig2 = openfig('NullSLERPResults/TopOpt/Gripper/FinalResults_NoOscillations/Monitoring_trust0d02_gJ5.fig');

iter1 = fig1.Children(end).Children.XData;
iter2 = fig2.Children(end).Children.XData;

cost1 = fig1.Children(end).Children.YData;
cost2 = fig2.Children(end).Children.YData;

constr1 = fig1.Children(end-2).Children.YData;
constr2 = fig2.Children(end-2).Children.YData;

lambda1 = fig1.Children(end-4).Children.YData;
lambda2 = fig2.Children(end-4).Children.YData;

eta1 = fig1.Children(end-7).Children.YData;
eta2 = fig2.Children(end-7).Children.YData;

close all;

colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]};

f1 = figure;
plot(iter1,cost1,iter2,cost2,'LineWidth',1.2);
grid on
grid minor
legend('$\eta^*=0.05$','$\eta^*=5.00$','interpreter','latex','FontSize', 12.5)
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f1.Position = [680   558   560/1.7   420];

f2 = figure;
plot(iter1,constr1,iter2,constr2,'LineWidth',1.2)
grid on
grid minor
legend('$\eta^*=0.05$','$\eta^*=5.00$','interpreter','latex','FontSize', 12.5)
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f2.Position = [680   558   560/1.7   420];

f3 = figure;
plot(iter1,lambda1,iter2,lambda2,'LineWidth',1.2)
grid on
grid minor
legend('$\eta^*=0.05$','$\eta^*=5.00$','interpreter','latex','FontSize', 12.5)
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f3.Position = [680   558   560/1.7   420];

f4 = figure;
h = plot(iter1(1:53),eta1(1:53),'--',iter1(53:end),eta1(53:end),iter2(1:48),eta2(1:48),'--',iter2(48:end),eta2(48:end),'LineWidth',1.2);
grid on
grid minor
legend('','$\eta^*=0.05$','','$\eta^*=5.00$','interpreter','latex','FontSize', 12.5)
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f4.Position = [680   558   560/1.7   420];
[h(1).Color, h(3).Color, h(5).Color] = colors{:};
[h(2).Color, h(4).Color, h(6).Color] = colors{:};



%% 1 Load case

clear;
close all;

fig1 = openfig('NullSLERPResults/TopOpt/MultiLoadBridge/FinalResults_NoOscillations/Monitoring_trust0d001_gJ10_1Loads.fig');

iter1 = fig1.Children(end).Children.XData;

cost1 = fig1.Children(end).Children.YData;

constr1 = fig1.Children(end-2).Children.YData;

lambda1 = fig1.Children(end-4).Children.YData;

eta1 = fig1.Children(end-7).Children.YData;

close all;

colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]};

f1 = figure;
plot(iter1,cost1,'LineWidth',1.2);
grid on
grid minor
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f1.Position = [680   558   560/1.7   420];

f2 = figure;
plot(iter1,constr1,'LineWidth',1.2)
grid on
grid minor
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f2.Position = [680   558   560/1.7   420];

f3 = figure;
plot(iter1,lambda1,'LineWidth',1.2)
grid on
grid minor
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f3.Position = [680   558   560/1.7   420];

f4 = figure;
h = plot(iter1(1:53),eta1(1:53),'--',iter1(54:end),eta1(54:end),'LineWidth',1.2);
grid on
grid minor
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f4.Position = [680   558   560/1.7   420];
h(1).Color = colors{1};
h(2).Color = colors{1};



%% 3 Loads case

clear;
close all;

fig1 = openfig('NullSLERPResults/TopOpt/MultiLoadBridge/FinalResults_NoOscillations/Monitoring_trust0d001_gJ10_3Loads.fig');

iter1 = fig1.Children(end).Children.XData;

cost1 = fig1.Children(end).Children.YData;

constr1 = fig1.Children(end-2).Children.YData;
constr2 = fig1.Children(end-3).Children.YData;
constr3 = fig1.Children(end-4).Children.YData;

lambda1 = fig1.Children(end-6).Children.YData;
lambda2 = fig1.Children(end-7).Children.YData;
lambda3 = fig1.Children(end-8).Children.YData;

eta1 = fig1.Children(end-11).Children.YData;

close all;

colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]};

f1 = figure;
plot(iter1,cost1,'LineWidth',1.2);
grid on
grid minor
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f1.Position = [680   558   560/1.7   420];
xlim([0 700])

f2 = figure;
plot(iter1,constr1,iter1,constr2,iter1,constr3,'LineWidth',1.2)
grid on
grid minor
legend('Load case 0','Load case 4','Load case 8','interpreter','latex','FontSize', 12.5)
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f2.Position = [680   558   560/1.7   420];
xlim([0 700])

f3 = figure;
plot(iter1,lambda1,iter1,lambda2,iter1,lambda3,'LineWidth',1.2)
grid on
grid minor
legend('Load case 0','Load case 4','Load case 8','interpreter','latex','FontSize', 12.5)
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f3.Position = [680   558   560/1.7   420];
xlim([0 700])

f4 = figure;
h = plot(iter1(1:53),eta1(1:53),'--',iter1(54:end),eta1(54:end),'LineWidth',1.2);
grid on
grid minor
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f4.Position = [680   558   560/1.7   420];
h(1).Color = colors{1};
h(2).Color = colors{1};
xlim([0 700])



%% All Loads case

clear;
close all;

fig1 = openfig('NullSLERPResults/TopOpt/MultiLoadBridge/FinalResults_NoOscillations/Monitoring_trust0d001_gJ10_9Loads.fig');

iter1 = fig1.Children(end).Children.XData;

cost1 = fig1.Children(end).Children.YData;

constr1 = fig1.Children(end-2).Children.YData;
constr2 = fig1.Children(end-3).Children.YData;
constr3 = fig1.Children(end-4).Children.YData;
constr4 = fig1.Children(end-5).Children.YData;
constr5 = fig1.Children(end-6).Children.YData;
constr6 = fig1.Children(end-7).Children.YData;
constr7 = fig1.Children(end-8).Children.YData;
constr8 = fig1.Children(end-9).Children.YData;
constr9 = fig1.Children(end-10).Children.YData;

lambda1 = fig1.Children(end-12).Children.YData;
lambda2 = fig1.Children(end-13).Children.YData;
lambda3 = fig1.Children(end-14).Children.YData;
lambda4 = fig1.Children(end-15).Children.YData;
lambda5 = fig1.Children(end-16).Children.YData;
lambda6 = fig1.Children(end-17).Children.YData;
lambda7 = fig1.Children(end-18).Children.YData;
lambda8 = fig1.Children(end-19).Children.YData;
lambda9 = fig1.Children(end-20).Children.YData;

eta1 = fig1.Children(end-23).Children.YData;

close all;

colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]};

f1 = figure;
plot(iter1,cost1,'LineWidth',1.2);
grid on
grid minor
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f1.Position = [680   558   560/1.7   420];

f2 = figure;
plot(iter1,constr1,iter1,constr2,iter1,constr3,iter1,constr4,iter1,constr5,iter1,constr6,iter1,constr7,iter1,constr8,'--',iter1,constr9,'--','LineWidth',1.2)
grid on
grid minor
legend('Load case 0','Load case 1','Load case 2','Load case 3','Load case 4','Load case 5','Load case 6','Load case 7','Load case 8','interpreter','latex','FontSize', 12.5)
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f2.Position = [680   558   560/1.7   420];

f3 = figure;
plot(iter1,lambda1,iter1,lambda2,iter1,lambda3,iter1,lambda4,iter1,lambda5,iter1,lambda6,iter1,lambda7,iter1,lambda8,'--',iter1,lambda9,'--','LineWidth',1.2)
grid on
grid minor
legend('Load case 0','Load case 1','Load case 2','Load case 3','Load case 4','Load case 5','Load case 6','Load case 7','Load case 8','interpreter','latex','FontSize', 12.5)
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f3.Position = [680   558   560/1.7   420];

f4 = figure;
h = plot(iter1(1:53),eta1(1:53),'--',iter1(54:end),eta1(54:end),'LineWidth',1.2);
grid on
grid minor
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f4.Position = [680   558   560/1.7   420];
h(1).Color = colors{1};
h(2).Color = colors{1};