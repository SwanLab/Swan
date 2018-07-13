function find_parameters

%% Problem circle
% Load data
load('data_circle.mat')

% Solver
eqn_circle = @(a,b,R) (aux_coords(:,1) - a).^2 + (aux_coords(:,2) - b).^2 - R^2;
fun = @(x) eqn_circle(x(1),x(2),x(3));

[x,fval,exitflag] = fsolve(fun,[1,1,1]);
fprintf('Center = [%g,%g];\n',x(1),x(2));
fprintf('Radius = %g;\n',x(3));

%% Problem line
% Load data
load('data_line.mat')

% Solver
eqn_line = @(m,n) aux_coords(:,2) - (m*aux_coords(:,1) + n);
fun = @(x) eqn_line(x(1),x(2));

[x,fval,exitflag] = fsolve(fun,[1,1,1]);
fprintf('m = %g;\n',x(1));
fprintf('n = %g;\n',x(2));

end