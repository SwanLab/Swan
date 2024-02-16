% CANTILEVER ACCELERATION CASE
% General config
s.problemType  = 'GENERAL';
s.geometryCase = 'CANTILEVER';
% Optimization config
s.tau     = 1e-1;
s.TOL     = 1e-3;
s.maxIter = 60;
s.lambda  = 10;
% Display final topology
s.result  = true;