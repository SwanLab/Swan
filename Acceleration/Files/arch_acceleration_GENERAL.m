% ARCH ACCELERATION CASE
% General config
s.problemType  = 'GENERAL';
s.geometryCase = 'ARCH';
% Optimization config
s.tau     = 5e-2;
s.TOL     = 1e-4;
s.maxIter = 60;
s.lambda  = 30;
% Display final topology
s.result  = true;