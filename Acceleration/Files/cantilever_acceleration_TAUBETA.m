% CANTILEVER ACCELERATION CASE
% General config
s.problemType  = 'TAU_BETA';
s.geometryCase = 'CANTILEVER';
% Optimization config
s.tau     = logspace(-3,-1,14);
s.beta    = 0:0.05:1;
s.TOL     = 1e-4;
s.maxIter = 80;
s.lambda  = 10;
s.gDescentType = 'Nesterov';
% Display final topology
s.result  = true;