sP.nLevels = 5;
sP.tol = 1e-8;
sP.maxIter = 1;
load('Por3D_RBmat.mat');
sP.nullSpace = RBbasisFree;

s.preconditioner = SmoothedAggregation(sP);
s.tol  = 1e-5;
s.type = 'PCG';