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

%[x1,res1]= run_pyamg_with_residuals2(A,b, 'B_nullspace', RBM);

s.type = 'ELASTIC';
s.LHS  = A;
x2     = pyAMG.create(s);