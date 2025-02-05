%% 3 gJs case

clear;
close all;

fig1 = openfig('NullSLERPResults/AcademicProblems/Delta0d02/Case3Eta0d2.fig');
fig2 = openfig('NullSLERPResults/AcademicProblems/Delta0d02/Case3Eta1.fig');
fig3 = openfig('NullSLERPResults/AcademicProblems/Delta0d02/Case3Eta5.fig');

iter1 = fig1.Children(end).Children.XData;
iter2 = fig2.Children(end).Children.XData;
iter3 = fig3.Children(end).Children.XData;

eta1 = fig1.Children(end-9).Children.YData;
eta2 = fig2.Children(end-9).Children.YData;
eta3 = fig3.Children(end-9).Children.YData;

close all;

colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]};

f1 = figure;
h = plot(iter1(1:1),eta1(1:1),'--',iter1(1:end),eta1(1:end),iter2(1:1),eta2(1:1),'--',iter2(1:end),eta2(1:end),iter3(1:1),eta3(1:1),'--',iter3(1:end),eta3(1:end),'LineWidth',1.2);
grid on
grid minor
legend('','$\eta^*=0.20$','','$\eta^*=1.00$','','$\eta^*=5.00$','interpreter','latex','FontSize', 12.5)
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f1.Position = [680   558   560/1.7   420];
[h(1).Color, h(3).Color, h(5).Color] = colors{:};
[h(2).Color, h(4).Color, h(6).Color] = colors{:};
%xlim([0 85])