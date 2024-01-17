% Gauss-Seidel Method in MATLAB
function [X, iter] = gauss_seidel(A, B)

% Separation of matrix A into lower triangular and upper triangular matrices
% A = D + L + U
D = diag(diag(A));
L = tril(A)- D;
U = triu(A)- D;

% allowable error in final answer
t = 10.0e-7;
iter = 0;
err = 1.0e10;

X = zeros(size(B));
while(norm(err) > t)
    iter = iter + 1;
    X_old = X;
    X = (D+L)\(B-U*X_old);% Gauss-Seidel formula
    err = X - X_old;% finding error
end
end