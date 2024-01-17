% Define the domain and grid spacing
Lx = 1;
Ly = 1;
nx = 31;
ny = 31;
hx = Lx / (nx - 1);
hy = Ly / (ny - 1);

% Define the exact solution
[X, Y] = meshgrid(0:hx:Lx, 0:hy:Ly);
u_exact = sin(pi*X) .* sin(pi*Y);

% Define the right-hand side function
f = -2*pi^2 * (sin(pi*X) .* sin(pi*Y));

% Define the boundary conditions
u = zeros(ny, nx);
u(1,:) = u_exact(1,:);  % Bottom boundary
u(ny,:) = u_exact(ny,:);  % Top boundary
u(:,1) = u_exact(:,1);  % Left boundary
u(:,nx) = u_exact(:,nx);  % Right boundary

% Solve the Poisson equation using Jacobi iteration
tol = 1e-6;
max_iter = 10000;
iter = 0;
residual = inf;
while residual > tol && iter < max_iter
    % Compute the residual
    r = zeros(ny, nx);
    for i = 2:ny-1
        for j = 2:nx-1
            r(i,j) = (u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) - 4*u(i,j)) / (hx^2 + hy^2) - f(i,j);
        end
    end
    residual = norm(r(:), 2);
   
    % Update the solution using Jacobi iteration
    for i = 2:ny-1
        for j = 2:nx-1
            u(i,j) = 0.25 * (u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) - hx^2*f(i,j));
        end
    end
   
    iter = iter + 1;
end

% Compute the error between the numerical and exact solutions
error = u - u_exact;

% Plot the numerical solution, exact solution, and error
subplot(1,3,1);
surf(X, Y, u);
xlabel('x');
ylabel('y');
zlabel('u');
title('Numerical solution');

subplot(1,3,2);
surf(X, Y, u_exact);
xlabel('x');
ylabel('y');
zlabel('u');
title('Exact solution');

subplot(1,3,3);
surf(X, Y, error);
xlabel('x');
ylabel('y');
zlabel('Error');
title('Error');