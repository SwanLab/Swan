function [exp1,exp2,exp3,exp4,exp5,exp6] = ExperimentingAcceleratedShapeOpt()
% 
exp1 = computeStandardCaseAdaptative();
exp2 = computeMomentumCase();
exp3 = computeConstantOneAdaptative();
% exp4 = computeStandardMMA();
% exp5 = computeMomentumCaseMMA();
% exp6 = computeConstantOneMMA();




figure()
hold on
plot(exp1.JV,'-+','LineWidth',2)
plot(exp2.JV,'-+','LineWidth',2)
plot(exp3.JV,'-+','LineWidth',2)
% plot(exp4.JV,'-+','LineWidth',2)
% plot(exp5.JV,'-+','LineWidth',2)
% plot(exp6.JV,'-+','LineWidth',2)
xlabel('Iteration','interpreter','latex')
ylabel('Cost','interpreter','latex')
box on
legend('PG - Standard','PG - Momentum','PG - Constant One','MMA - Standard','MMA - Momentum','MMA - Constant One','interpreter','latex')
set(gca,'FontSize',14,'TickLabelInterpreter','latex')
hold off
figure()
hold on
plot(exp1.incXvalues,'-+','LineWidth',2)
plot(exp2.incXvalues,'-+','LineWidth',2)
plot(exp3.incXvalues,'-+','LineWidth',2)
% plot(exp4.incXvalues,'-+','LineWidth',2)
% plot(exp5.incXvalues,'-+','LineWidth',2)
% plot(exp6.incXvalues,'-+','LineWidth',2)
xlabel('Iteration','interpreter','latex')
ylabel('incX','interpreter','latex')
box on
legend('PG - Standard','PG - Momentum','PG - Constant One','MMA - Standard','MMA - Momentum','MMA - Constant One','interpreter','latex')
set(gca,'FontSize',14,'TickLabelInterpreter','latex')
hold off
figure()
hold on
plot(exp1.betaV,'-+','LineWidth',2)
plot(exp2.betaV,'-+','LineWidth',2)
plot(exp3.betaV,'-+','LineWidth',2)
xlabel('Iteration','interpreter','latex')
ylabel('$\beta_k$','interpreter','latex')
box on
legend('Standard','Momentum','Constant One','interpreter','latex')
set(gca,'FontSize',14,'TickLabelInterpreter','latex')
hold off
figure()
hold on
plot(exp1.tV,'-+','LineWidth',2)
plot(exp2.tV,'-+','LineWidth',2)
plot(exp3.tV,'-+','LineWidth',2)
xlabel('Iteration','interpreter','latex')
ylabel('Line search ($\tau$)','interpreter','latex')
box on
legend('PG - Standard','PG - Momentum','PG - Constant One','interpreter','latex')
set(gca,'FontSize',14,'TickLabelInterpreter','latex')
end

function plotSubFigure(nFigures,nF,name,y1,y2,y3)
subplot(2,2,nF)
plot(y1,'-+','LineWidth',3);
hold on
plot(y2,'-+','LineWidth',3);
hold on
plot(y3,'-+','LineWidth',3);
legend('Standard','Momentum','Constant One')
title(name)
end


function solver = computeMomentumCase()
s.momentumParams.type = 'CONVEX';
s.maxIter = 100;
s.TOL = 1e-12;
solver = ShapeOptimizationSolver(s);
solver.solve();
end

function solver = computeMomentumCaseAdaptative()
s.momentumParams.type = 'CONVEX';
s.maxIter = 100;
s.TOL = 1e-12;
solver = ShapeOptimizationSolver(s);
solver.solveAdaptative();
end


function solver = computeMomentumCaseMMA()
s.momentumParams.type = 'CONVEX';
s.maxIter = 100;
s.TOL = 1e-12;
solver = MMAShapeOptimizationSolver(s);
solver.solve();
end

function solver = computeStandardCase()
s.momentumParams.type = 'CONSTANT';
s.momentumParams.value = 0;
s.maxIter = 100;
s.TOL = 1e-12;
solver = ShapeOptimizationSolver(s);
solver.solve();
end

function solver = computeStandardCaseAdaptative()
s.momentumParams.type = 'CONSTANT';
s.momentumParams.value = 0;
s.maxIter = 100;
s.TOL = 1e-12;
solver = ShapeOptimizationSolver(s);
solver.solveAdaptative();
end

function solver = computeStandardMMA()
s.momentumParams.type = 'CONSTANT';
s.momentumParams.value = 0;
s.maxIter = 100;
s.TOL = 1e-12;
solver = MMAShapeOptimizationSolver(s);
solver.solve();
end

function solver = computeConstantOne()
s.momentumParams.type = 'CONSTANT';
s.momentumParams.value = 1;
s.maxIter = 100;
s.TOL = 1e-12;
solver = ShapeOptimizationSolver(s);
solver.solve();
end

function solver = computeConstantOneAdaptative()
s.momentumParams.type = 'CONSTANT';
s.momentumParams.value = 0.3;
s.maxIter = 100;
s.TOL = 1e-12;
solver = ShapeOptimizationSolver(s);
solver.solveAdaptative();
end

function solver = computeConstantOneMMA()
s.momentumParams.type = 'CONSTANT';
s.momentumParams.value = 1;
s.maxIter = 100;
s.TOL = 1e-12;
solver = MMAShapeOptimizationSolver(s);
solver.solve();
end
