function ExperimentingAcceleratedShapeOpt()




exp1 = computeStandardCase();
exp2 = computeMomentumCase();
exp3 = computeConstantOne();


nFigures = 5;
fh = figure('units', 'pixels');
hold on
fh.set('Position',[4000 1500 3000 500])


plotSubFigure(nFigures,1,'Cost',exp1.JV,exp2.JV,exp3.JV)
plotSubFigure(nFigures,2,'LineSearch',exp1.tV,exp2.tV,exp3.tV)
plotSubFigure(nFigures,3,'Beta',exp1.betaV,exp2.betaV,exp3.betaV)
plotSubFigure(nFigures,4,'IncX',exp1.incXvalues,exp2.incXvalues,exp3.incXvalues)


end

function plotSubFigure(nFigures,nF,name,y1,y2,y3)
subplot(1,nFigures,nF)
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

function solver = computeStandardCase()
s.momentumParams.type = 'CONSTANT';
s.momentumParams.value = 0;
s.maxIter = 100;
s.TOL = 1e-12;
solver = SimpleShapeOptimizationSolver(s);
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
