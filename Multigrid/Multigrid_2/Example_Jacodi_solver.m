% Define grid and boundary conditions
nx = 100;
ny = 100;
Lx = 1.0;
Ly = 1.0;
dx = Lx / (nx-1);
dy = Ly / (ny-1);
x = linspace(0, Lx, nx);
y = linspace(0, Ly, ny);
[X, Y] = meshgrid(x, y);

% Set boundary conditions
u = zeros(nx, ny);
u(:,1) = 0;    % Dirichlet boundary condition on left edge
u(:,end) = 1;  % Dirichlet boundary condition on right edge

% Set up coefficient matrix and right-hand side vector
A = gallery('poisson', nx, ny);
b = reshape(A*u(:), [], 1);

% Apply Dirichlet boundary conditions to matrix and vector
ind = [1:nx, (ny-1)*nx+1:ny*nx];   % indices of boundary points
A(ind,:) = sparse(1:numel(ind), ind, ones(numel(ind),1), numel(ind), nx*ny);
b(ind) = u(ind);

% Jacobi smoothing iteration
tol = 1e-7;
err = tol+1;
iter = 0;
while err > tol
    unew = u;
    for i = 2:nx-1
        for j = 2:ny-1
            unew(i,j) = 0.25*(u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) - dx^2*b((i-1)*ny+j));
        end
    end
    err = norm(unew - u, 'inf') / norm(unew, 'inf');
    u = unew;
    iter = iter + 1;
end

% Plot solution
figure;
surf(X, Y, u');
xlabel('x');
ylabel('y');
zlabel('u');
title(['Poisson equation solution (Jacobi), iter = ' num2str(iter)]);
