%% Simple topology optimization example
% Note that the beta term

s.momentumParams.type = 'CONSTANT';
s.momentumParams.value = 0;
s.maxIter = 100;
s.TOL = 1e-12;
solver = SimpleShapeOptimizationSolver(s);
solver.solve();