%% 3 gJs case

clear;
close all;

fig1 = openfig('NullSLERPResults/AcademicProblems/Delta0d02/Case4Eta0d2.fig');
fig2 = openfig('NullSLERPResults/AcademicProblems/Delta0d02/Case4Eta1.fig');
fig3 = openfig('NullSLERPResults/AcademicProblems/Delta0d02/Case4Eta5.fig');

iter1 = fig1.Children(end).Children.XData;
iter2 = fig2.Children(end).Children.XData;
iter3 = fig3.Children(end).Children.XData;

eta1 = fig1.Children(end-9).Children.YData;
eta2 = fig2.Children(end-9).Children.YData;
eta3 = fig3.Children(end-9).Children.YData;

etaMax1 = 1./fig1.Children(end-7).Children.YData;
etaMax2 = 1./fig2.Children(end-7).Children.YData;
etaMax3 = 1./fig3.Children(end-7).Children.YData;

close all;

colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]};

f1 = figure;
h = plot(iter1,eta1,iter1,etaMax1,'--',iter2,eta2,iter2,etaMax2,'--',iter3,eta3,iter3,etaMax3,'--','LineWidth',1.2);
grid on
grid minor
legend('$\eta^*=0.20$','$1/\tau_{k(0.20)}^*$','$\eta^*=1.00$','$1/\tau_{k(1.00)}^*$','$\eta^*=5.00$','$1/\tau_{k(5.00)}^*$','interpreter','latex','FontSize', 12.5)
xlabel('Iteration','Interpreter','latex','FontSize', 12.5)
f1.Position = [680   558   560/1.7   420];
[h(1).Color, h(3).Color, h(5).Color] = colors{:};
[h(2).Color, h(4).Color, h(6).Color] = colors{:};
f1.Children(2).YScale = 'linear';
%xlim([0 180])
%ylim([0 600])