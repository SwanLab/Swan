% Material parameters
mu = 1;
lambda = 1;

% Matrices

u = [1 1 1];

F = eye(3) + Grad(u); % deformation gradient
jac = det(F);
cofac = adjoint(F)';
% invFt = inv(F);
invFt = cofac/jac;

piola = mu*(F-invFt) + lambda*(jac-1)*jac*invFt;