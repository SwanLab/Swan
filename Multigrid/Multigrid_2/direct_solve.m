function u = direct_solve(f, h)
% Solves Poisson equation directly on a square domain with Dirichlet BCs using a 5-point stencil
% Input: f - right-hand side function (L-1 by L-1 matrix)
%        h - grid spacing
% Output: u - numerical solution (L-1 by L-1 matrix)

% Define matrix of coefficients
L = size(f, 1);
A = sparse(L^2, L^2);
for j = 1:L
    for i = 1:L
        row = (j-1)*L + i;
        if i == 1 || i == L || j == 1 || j == L
            A(row, row) = 1;
        else
            A(row, row) = -4;
            A(row, row-1) = 1;
            A(row, row+1) = 1;
            A(row, row-L) = 1;
            A(row, row+L) = 1;
        end
    end
end

% Convert matrix to full format and solve
A = full(A);
f = reshape(f', [], 1); % reshape f into a column vector
u = A \ (h^2 * f); % solve linear system
u = reshape(u, L, [])'; % reshape u into a matrix
disp(u)
end
