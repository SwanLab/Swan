% CANTILEVER ACCELERATION CASE
% General config
s.problemType  = 'SOLVER';
s.geometryCase = 'CANTILEVER3D';
s.geometryType = '3D';
s.matrixFree = false;
% Optimization config
s.tau     = 1e-2;
s.solverType = 'CONJUGATE GRADIENT';%'DIRECT';%
s.maxTol = 1e-1;
s.minTol = 1e-3;
s.eta    = 2;
s.TOL     = 5e-3;
s.maxIter = 5e3;
s.lambda  = 60;
% Display final topology
s.result  = true;