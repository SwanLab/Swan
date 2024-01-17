

function A = poisson2d(L)
% Creates the 2D Poisson matrix with Dirichlet boundary conditions
% on an n x n grid.
%
% Inputs:
% n - The number of interior grid points in each dimension (i.e. not including
%     boundary points)
%
% Outputs:
% A - The n^2 x n^2 Poisson matrix

% Number of interior grid points in each dimension
m = L-1;

% Define the matrix A
A = zeros(m^2);

% Fill in the interior grid points
for i = 2:m
    for j = 2:m
        index = (i-1)*m + j;
        A(index, index) = -4;
        A(index, index-1) = 1;
        A(index, index+1) = 1;
        A(index, index-m) = 1;
        A(index, index+m) = 1;
    end
end

% Fill in the boundary conditions (Dirichlet conditions)
for i = 1:m
    index = i;
    A(index, :) = 0;
    A(index, index) = 1;
    
    index = m*(m-1) + i;
    A(index, :) = 0;
    A(index, index) = 1;
    
    index = (i-1)*m + 1;
    A(index, :) = 0;
    A(index, index) = 1;
    
    index = i*m;
    A(index, :) = 0;
    A(index, index) = 1;
end

end

