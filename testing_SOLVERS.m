% TO TEST SOLVERS
close all
filename = 'cantilever_acceleration_SOLVER';
filename2 = 'cantilever_acceleration_SOLVER_2';

a = AccelerationExperiments(filename);
a2 = AccelerationExperiments(filename2);

p1 = a.experiment.problem;
p2 = a2.experiment.problem;
figure()
hold on
plot(p1.J,'-+','LineWidth',1.5)
plot(p2.J,'-+','LineWidth',1.5)
xlabel('Iteration','interpreter','latex')
ylabel('$J$','interpreter','latex')
legend('Direct','CG + $\beta_{adapt.}$','interpreter','latex')
set(gca,'FontSize',14,'TickLabelInterpreter','latex')
set(gca, 'YScale', 'log')
xlim([1, max(numel(p1.J),numel(p2.J))])
grid minor
box on
hold off

figure()
hold on
plot(p1.costFields(1,:),'-+','LineWidth',1.5)
plot(p2.costFields(1,:),'-+','LineWidth',1.5)
xlabel('Iteration','interpreter','latex')
ylabel('Compliance','interpreter','latex')
set(gca,'FontSize',14,'TickLabelInterpreter','latex')
legend('Direct','CG + $\beta_{adapt.}$','interpreter','latex')
set(gca, 'YScale', 'log')
xlim([1, max(numel(p1.J),numel(p2.J))])
grid minor
box on
hold off

figure()
hold on
plot(p1.costFields(2,:),'-+','LineWidth',1.5)
plot(p2.costFields(2,:),'-+','LineWidth',1.5)
xlabel('Iteration','interpreter','latex')
ylabel('Volume frac.','interpreter','latex')
set(gca,'FontSize',14,'TickLabelInterpreter','latex')
legend('Direct','CG + $\beta_{adapt.}$','interpreter','latex')
set(gca, 'YScale', 'log')
xlim([1, max(numel(p1.J),numel(p2.J))])
grid minor
box on
hold off
