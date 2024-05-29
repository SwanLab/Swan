% CANTILEVER ACCELERATION CASE
% General config
s.problemType  = 'GENERAL';
s.geometryCase = 'CANTILEVER';
% Optimization config
s.tau     = 1e-1;
s.TOL     = 5e-4;
s.maxIter = 80;
s.lambda  = 10;
% Display final topology
s.result  = true;