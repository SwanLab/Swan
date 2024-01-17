

function u_exact = exact_soultion(x,y)

% Compute the exact solution for the Poisson equation
% for f(x,y) one;
% We still assume Homogeneous Dirichlet BCs on a 2D grid

    u_exact = sin(pi*x).*sin(pi*y);

end