% BRIDGE ACCELERATION CASE
% General config
s.problemType  = 'GENERAL';
s.geometryCase = 'BRIDGE';
% Optimization config
s.tau     = 5e-2;
s.TOL     = 1e-4;
s.maxIter = 60;
s.lambda  = 10;
% Display final topology
s.result  = true;