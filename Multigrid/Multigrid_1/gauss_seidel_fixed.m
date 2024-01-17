% Gauss-Seidel Method in MATLAB
function [X, conv, local_iter] = gauss_seidel_fixed(A ,B, X0, num_gs_steps)

% Separation of matrix A into lower triangular and upper triangular matrices
% A = D + L + U
D = diag(diag(A));
L = tril(A)- D;
U = triu(A)- D;

% allowable error in final answer
t = 10.0e-7;

X = X0;
conv = 0;
local_iter = 0;
for i =1:num_gs_steps
    X_old = X;
    X = (D+L)\(B-U*X_old);% Gauss-Seidel formula
    err = X - X_old;% finding error
    local_iter = local_iter + 1;
    if(norm(err)<t)
        conv = 1;
        return
    end
end