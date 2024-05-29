function [exp1,exp2,exp3,exp4,exp5] = ExperimentingAcceleratedShapeOpt()
% 
exp1 = computeStandardCase();
exp2 = computeMomentumCase();
exp3 = computeConstantOne();

exp4 = computeNesterovMomentumCase();
exp5 = computeNesterovConstantOne();
exp4 = computeStandardMMA();
% exp5 = computeMomentumCaseMMA();
% exp6 = computeConstantOneMMA();




figure()
hold on
plot(exp1.JV,'-+','LineWidth',2)
plot(exp2.JV,'-+','LineWidth',2)
plot(exp3.JV,'-+','LineWidth',2)
plot(exp4.JV,'-+','LineWidth',2)
plot(exp5.JV,'-+','LineWidth',2)
% plot(exp6.JV,'-+','LineWidth',2)
xlabel('Iteration','interpreter','latex')
ylabel('Cost','interpreter','latex')
set(gca, 'YScale', 'log')
box on
% legend('PG - Standard','PG - Momentum','PG - Constant One','MMA - Standard','MMA - Momentum','MMA - Constant One','interpreter','latex')
legend('PG - Standard','Momentum Poly.','Constant One Poly.','Momentum Nest.','Constant one Nest.','interpreter','latex')
set(gca,'FontSize',14,'TickLabelInterpreter','latex')
grid minor
hold off
figure()
hold on
plot(exp1.incJV,'-+','LineWidth',2)
plot(exp2.incJV,'-+','LineWidth',2)
plot(exp3.incJV,'-+','LineWidth',2)
plot(exp4.incJV,'-+','LineWidth',2)
plot(exp5.incJV,'-+','LineWidth',2)
% plot(exp6.incJV,'-+','LineWidth',2)
xlabel('Iteration','interpreter','latex')
ylabel('Cost increment ($\Delta J$)','interpreter','latex')
set(gca, 'YScale', 'log')
box on
% legend('PG - Standard','PG - Momentum','PG - Constant One','MMA - Standard','MMA - Momentum','MMA - Constant One','interpreter','latex')
legend('PG - Standard','Momentum Poly.','Constant One Poly.','Momentum Nest.','Constant one Nest.','interpreter','latex')
set(gca,'FontSize',14,'TickLabelInterpreter','latex')
grid minor
hold off

% --- 
minJV = [min(exp1.JV);min(exp2.JV);min(exp3.JV);min(exp4.JV);min(exp5.JV)];
figure()
hold on
plot(exp1.JV-minJV,'-+','LineWidth',2)
plot(exp2.JV-minJV,'-+','LineWidth',2)
plot(exp3.JV-minJV,'-+','LineWidth',2)
plot(exp4.JV-minJV,'-+','LineWidth',2)
plot(exp5.JV-minJV,'-+','LineWidth',2)
% plot(exp6.JV,'-+','LineWidth',2)
xlabel('Iteration','interpreter','latex')
ylabel('Cost','interpreter','latex')
set(gca, 'YScale', 'log')
box on
% legend('PG - Standard','PG - Momentum','PG - Constant One','MMA - Standard','MMA - Momentum','MMA - Constant One','interpreter','latex')
legend('PG - Standard','Momentum Poly.','Constant One Poly.','Momentum Nest.','Constant one Nest.','interpreter','latex')
set(gca,'FontSize',14,'TickLabelInterpreter','latex')
grid minor
hold off
% figure()
% hold on
% plot(exp1.betaV,'-+','LineWidth',2)
% plot(exp2.betaV,'-+','LineWidth',2)
% plot(exp3.betaV,'-+','LineWidth',2)
% xlabel('Iteration','interpreter','latex')
% ylabel('$\beta_k$','interpreter','latex')
% box on
% legend('Standard','Momentum','Constant One','interpreter','latex')
% set(gca,'FontSize',14,'TickLabelInterpreter','latex')
% hold off
% figure()
% hold on
% plot(exp1.tV,'-+','LineWidth',2)
% plot(exp2.tV,'-+','LineWidth',2)
% plot(exp3.tV,'-+','LineWidth',2)
% xlabel('Iteration','interpreter','latex')
% ylabel('Line search ($\tau$)','interpreter','latex')
% box on
% legend('PG - Standard','PG - Momentum','PG - Constant One','interpreter','latex')
% set(gca,'FontSize',14,'TickLabelInterpreter','latex')
end

% function plotSubFigure(nFigures,nF,name,y1,y2,y3)
% subplot(2,2,nF)
% plot(y1,'-+','LineWidth',3);
% hold on
% plot(y2,'-+','LineWidth',3);
% hold on
% plot(y3,'-+','LineWidth',3);
% legend('Standard','Momentum','Constant One')
% title(name)
% end


function solver = computeNesterovMomentumCase()
s.momentumParams.type = 'CONVEX';
s.maxIter = 100;
s.TOL = 5e-4;
s.tau = 80;
solver = ShapeOptimizationSolver(s);
solver.solveNesterov();
end

function solver = computeMomentumCase()
s.momentumParams.type = 'CONVEX';
s.maxIter = 100;
s.TOL = 5e-4;
s.tau = 80;
solver = ShapeOptimizationSolver(s);
solver.solve();
end

function solver = computeMomentumCaseAdaptative()
s.momentumParams.type = 'CONVEX';
s.maxIter = 100;
s.TOL = 5e-4;
s.tau = 80;
solver = ShapeOptimizationSolver(s);
solver.solveAdaptative();
end


function solver = computeMomentumCaseMMA()
s.momentumParams.type = 'CONVEX';
s.maxIter = 100;
s.TOL = 5e-4;
solver = MMAShapeOptimizationSolver(s);
solver.solve();
end

function solver = computeStandardCase()
s.momentumParams.type = 'CONSTANT';
s.momentumParams.value = 0;
s.maxIter = 100;
s.TOL = 5e-4;
s.tau = 100;
solver = ShapeOptimizationSolver(s);
solver.solve();
end

function solver = computeStandardCaseAdaptative()
s.momentumParams.type = 'CONSTANT';
s.momentumParams.value = 0;
s.maxIter = 100;
s.TOL = 5e-4;
s.tau = 80;
solver = ShapeOptimizationSolver(s);
solver.solveAdaptative();
end

function solver = computeStandardMMA()
s.momentumParams.type = 'CONSTANT';
s.momentumParams.value = 0;
s.maxIter = 100;
s.TOL = 5e-4;
solver = MMAShapeOptimizationSolver(s);
solver.solve();
end

function solver = computeConstantOne()
s.momentumParams.type = 'CONSTANT';
s.momentumParams.value = 1;
s.maxIter = 100;
s.TOL = 5e-4;
s.tau = 80;
solver = ShapeOptimizationSolver(s);
solver.solve();
end

function solver = computeNesterovConstantOne()
s.momentumParams.type = 'CONSTANT';
s.momentumParams.value = 1;
s.maxIter = 100;
s.TOL = 5e-4;
s.tau = 80;
solver = ShapeOptimizationSolver(s);
solver.solveNesterov();
end

function solver = computeConstantOneAdaptative()
s.momentumParams.type = 'CONSTANT';
s.momentumParams.value = 0.3;
s.maxIter = 100;
s.TOL = 5e-4;
s.tau = 80;
solver = ShapeOptimizationSolver(s);
solver.solveAdaptative();
end

function solver = computeConstantOneMMA()
s.momentumParams.type = 'CONSTANT';
s.momentumParams.value = 1;
s.maxIter = 100;
s.TOL = 5e-4;
solver = MMAShapeOptimizationSolver(s);
solver.solve();
end
