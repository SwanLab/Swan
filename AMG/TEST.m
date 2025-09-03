clear all
load('Por3D_LHS.mat');
A=LHS;
load('Por3D_RHS.mat');
b=RHS;
load('Por3D_RBmat.mat');
if exist("RBbasisFree")
    RBM = RBbasisFree;
else
    RBM = B;
end

[x1,res1]= run_pyamg_with_residuals2(A,b, 'B_nullspace', RBM);

s.type = 'ELASTIC';
s.LHS  = A;
s.nullSpace = RBM;
s.nLevels = 5;
s.tol = 1e-8;
s.maxIter = 1;
p     = pyAMG.create(s);
%x2 = p.solve(b);

P = @(r) p.solve(r);
solver = PCG();
[x3,res3] = solver.solve(@(x) A*x,b,zeros(size(b)),P,1e-8);