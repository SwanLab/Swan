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
% RBM = B;
% u = elasticity_amg_solver(A,B,RBM);
[x,res]= run_pyamg_with_residuals2(A,b, 'B_nullspace', RBM);