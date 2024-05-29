function computeAdaptativeLineSearch()
    exp1 = computeStandardCase();
    exp2 = computeMomentumCase();
    exp3 = computeStandardCaseAdaptative();
    exp4 = computeMomentumCaseAdaptative();
    exp5 = computeConstantOneAdaptative();
    figure()
    hold on
    plot(exp1.JV,'-+','LineWidth',2)
    plot(exp2.JV,'-+','LineWidth',2)
    plot(exp3.JV,'-+','LineWidth',2)
    plot(exp4.JV,'-+','LineWidth',2)
    plot(exp5.JV,'-+','LineWidth',2)
    xlabel('Iteration','interpreter','latex')
    ylabel('Cost','interpreter','latex')
    box on
    legend('PG - Standard','PG - Momentum','PG - Standard adaptative','PG - Momentum adaptative','PG - ct. 1 adaptative','interpreter','latex')
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
    xlabel('Iteration','interpreter','latex')
    ylabel('Cost increment ($\Delta J$)','interpreter','latex')
    set(gca, 'YScale', 'log')
    box on
    legend('PG - Standard','PG - Momentum','PG - Standard adaptative','PG - Momentum adaptative','PG - ct. 1 adaptative','interpreter','latex')
    set(gca,'FontSize',14,'TickLabelInterpreter','latex')
    grid minor
    hold off
    figure()
    hold on
    plot(exp1.tV,'-+','LineWidth',2)
    plot(exp2.tV,'-+','LineWidth',2)
    plot(exp3.tV,'-+','LineWidth',2)
    plot(exp4.tV,'-+','LineWidth',2)
    plot(exp5.tV,'-+','LineWidth',2)
    xlabel('Iteration','interpreter','latex')
    ylabel('Line-search ($\tau$)','interpreter','latex')
    box on
    legend('PG - Standard','PG - Momentum','PG - Standard adaptative','PG - Momentum adaptative','PG - ct. 1 adaptative','interpreter','latex')
    set(gca,'FontSize',14,'TickLabelInterpreter','latex')
    set(gca, 'YScale', 'log')
    grid minor
    hold off
end

function solver = computeStandardCase()
    s.momentumParams.type = 'CONSTANT';
    s.momentumParams.value = 0;
    s.maxIter = 300;
    s.TOL = 1e-3;
    s.tau = 100;
    solver = ShapeOptimizationSolver(s);
    solver.solve();
end

function solver = computeStandardCaseAdaptative()
    s.momentumParams.type = 'CONSTANT';
    s.momentumParams.value = 0;
    s.maxIter = 100;
    s.TOL = 1e-3;
    s.tau = 10;
    solver = ShapeOptimizationSolver(s);
    solver.solveAdaptative();
end

function solver = computeConstantOneAdaptative()
    s.momentumParams.type = 'CONSTANT';
    s.momentumParams.value = 0.3;
    s.maxIter = 100;
    s.TOL = 1e-3;
    s.tau = 10;
    solver = ShapeOptimizationSolver(s);
    solver.solveAdaptative();
end

function solver = computeMomentumCaseAdaptative()
    s.momentumParams.type = 'CONVEX';
    s.maxIter = 100;
    s.TOL = 1e-3;
    s.tau = 10;
    solver = ShapeOptimizationSolver(s);
    solver.solveAdaptative();
end

function solver = computeMomentumCase()
s.momentumParams.type = 'CONVEX';
s.maxIter = 100;
s.TOL = 1e-3;
s.tau = 100;
solver = ShapeOptimizationSolver(s);
solver.solve();
end
