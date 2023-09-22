function [exp1,exp2,exp3,exp4,exp5,exp6] = ExperimentingAcceleratedShapeOpt()
% 
exp1 = computeStandardCase();
exp2 = computeMomentumCase();
exp3 = computeConstantOne();
exp4 = computeStandardMMA();
exp5 = computeMomentumCaseMMA();
exp6 = computeConstantOneMMA();

% nFigures = 3;
% % fh = figure('units', 'pixels');
% figure()
% hold on
% fh.set('Position',[4000 1500 3000 500])


% plotSubFigure(nFigures,1,'Cost',exp1.JV,exp2.JV,exp3.JV)
% plotSubFigure(nFigures,2,'LineSearch',exp1.tV,exp2.tV,exp3.tV)
% plotSubFigure(nFigures,3,'Beta',exp1.betaV,exp2.betaV,exp3.betaV)
% plotSubFigure(nFigures,4,'IncX',exp1.incXvalues,exp2.incXvalues,exp3.incXvalues)

% plotSubFigure(nFigures,1,'Cost',exp4.JV,exp5.JV,exp6.JV)
% plotSubFigure(nFigures,2,'Beta',exp4.betaV,exp5.betaV,exp6.betaV)
% plotSubFigure(nFigures,3,'IncX',exp4.incXvalues,exp5.incXvalues,exp6.incXvalues)


figure()
hold on
plot(exp1.JV,'-+','LineWidth',2)
plot(exp2.JV,'-+','LineWidth',2)
plot(exp3.JV,'-+','LineWidth',2)
plot(exp4.JV,'-+','LineWidth',2)
plot(exp5.JV,'-+','LineWidth',2)
plot(exp6.JV,'-+','LineWidth',2)
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
plot(exp4.incXvalues,'-+','LineWidth',2)
plot(exp5.incXvalues,'-+','LineWidth',2)
plot(exp6.incXvalues,'-+','LineWidth',2)
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

function solver = computeConstantOneMMA()
s.momentumParams.type = 'CONSTANT';
s.momentumParams.value = 1;
s.maxIter = 100;
s.TOL = 1e-12;
solver = MMAShapeOptimizationSolver(s);
solver.solve();
end
