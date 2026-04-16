
clear;
close all;

fig = openfig('MonitoringDensityVerticalCantileverGlobalCircle0.6.fig');

iter = fig.Children(end).Children.XData(1:301);

compl = fig.Children(end-1).Children.YData(1:301);
per = fig.Children(end-2).Children.YData(1:301);

volConstr = fig.Children(end-3).Children.YData(1:301);
perConstr = fig.Children(end-4).Children.YData(1:301);

lPer = fig.Children(end-7).Children.YData(1:301);

close all;

colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]};
cmap = colororder();

f1 = figure;
plot(iter,compl,iter,per,'LineWidth',1.2);
grid on
grid minor
legend('Compliance','Neg IP penalty','interpreter','latex','FontSize', 12.5)
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f1.Position = [680   558   560/1.7   420];

f2 = figure;
plot(iter,volConstr,iter,perConstr,'LineWidth',1.2)
grid on
grid minor
legend('Volume','Perimeter','interpreter','latex','FontSize', 12.5)
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f2.Position = [680   558   560/1.7   420];

f3 = figure;
plot(iter,lPer,'color', cmap(2,:),'LineWidth',1.2)
grid on
grid minor
legend('Perimeter','interpreter','latex','FontSize', 12.5)
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f3.Position = [680   558   560/1.7   420];