clear
clc

% =========================================================================
% 1) CREATE MESH
% =========================================================================
mesh = UnitTriangleMesh(100,100);

% =========================================================================
% 2) CREATE FINITE ELEMENT FUNCTIONS (velocity u, pressure p)
% =========================================================================
u = LagrangianFunction.create(mesh, 2, 'P2');
p = LagrangianFunction.create(mesh, 1, 'P1');
n_dofs = u.nDofs + p.nDofs;

% =========================================================================
% 3) SET UP BOUNDARY CONDITIONS (Dirichlet)
% =========================================================================
isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1))) < 1e-12);
isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1))) < 1e-12);
isBottom = @(coor) (abs(coor(:,2) - min(coor(:,2))) < 1e-12);
isTop    = @(coor) (abs(coor(:,2) - max(coor(:,2))) < 1e-12);

dir_vel{1}.domain    = @(coor) isTop(coor);
dir_vel{1}.direction = [1,2];
dir_vel{1}.value     = [100,0];

dir_vel{2}.domain    = @(coor) isBottom(coor) | isRight(coor) | isLeft(coor);
dir_vel{2}.direction = [1,2];
dir_vel{2}.value     = [0,0]; 

dirichlet = [];
dir_dofs = [];
for i = 1:length(dir_vel)
    dirDofs = u.getDofsFromCondition(dir_vel{i}.domain);
    nodes = 1 + (dirDofs(2:2:end) - 2) / u.ndimf;
    iNod = repelem(nodes, 2)';
    dofType = repmat((1:2)', length(nodes), 1);
    valmat = repmat(dir_vel{i}.value', length(nodes), 1);
    dirichlet = [dirichlet; [iNod', dofType, valmat]];
    dir_dofs = [dir_dofs; dirDofs(:)];
end

% =========================================================================
% 4) APPLIED FORCES
% =========================================================================
sAF.fHandle = @(coor) [0.*coor, 0.*coor];
sAF.ndimf   = 2;
sAF.mesh    = mesh;
forcesFormula = AnalyticalFunction(sAF);

% =========================================================================
% 5) BUILD MATRICES FOR THE LINEAR STOKES PROBLEM
% =========================================================================
% 1) Diffusion (viscous) part
s1.type  = 'StiffnessMatrix';
s1.mesh  = mesh;
s1.test  = u;
s1.trial = u;
LHS1 = LHSIntegrator.create(s1);
K = LHS1.compute();

% 2) Divergence part for incompressibility
s2.type = 'WeakDivergence';
s2.mesh = mesh;
s2.trial = p;
s2.test  = u;
LHS2 = LHSIntegrator.create(s2);
D = LHS2.compute();

% 3) Nonlinear convective term
s3.type = 'NonLinearNS';
s3.mesh = mesh;
s3.trial = u;
s3.test = u;
LHS3 = LHSIntegrator.create(s3);

% 3) Forces vector
s4.type          = 'Stokes';
s4.mesh          = mesh;
s4.velocityFun   = u;
s4.pressureFun   = p;
s4.forcesFormula = forcesFormula;
RHSint = RHSIntegrator.create(s4);
F = RHSint.integrate();

sz = size(D, 2);
BB = sparse(sz, sz);

uD = dirichlet(:,3);
LHS = [K, D; D', BB];

pFixDof  = u.nDofs+1;
LHS(pFixDof, :) = 0;
LHS(:, pFixDof) = 0;
LHS(pFixDof, pFixDof) = 1;

R  = -LHS*sparse(dir_dofs,1,uD,n_dofs,1);
RHS = F + R;
RHS(pFixDof) = 0;

% =========================================================================
% 6) SOLVE SYSTEM
% =========================================================================
maxIter = 100;        % maximum Picard iterations
tol     = 1e-8;      % tolerance for convergence
x       = zeros(n_dofs,1); % initial guess
lambda  = 1;       % optional relaxation factor for Picard
free_dofs = setdiff(1:n_dofs, dir_dofs);

for iter = 1:maxIter

    C = LHS3.compute(u);
    
    % Build the linear system for this iteration
    LHS = [K + C, D; 
           D',  BB];
    LHS(pFixDof, :) = 0;
    LHS(:, pFixDof) = 0;
    LHS(pFixDof, pFixDof) = 1;   

    % Rebuild boundary offset (because LHS changed)
    R = -LHS(:,dir_dofs) * uD;
    RHS = F + R - LHS * x;
    RHS(pFixDof) = 0;
    
    % Restrict to free dofs and solve
    LHSr = LHS(free_dofs, free_dofs);
    RHSr = RHS(free_dofs);
    dx = LHSr\RHSr;
    
    % Update free-dof solution with relaxation
    x(free_dofs) = x(free_dofs) + lambda*dx;
    
    % Enforce Dirichlet dofs
    fullx = x;
    if ~isempty(dir_dofs)
        fullx(dir_dofs) = uD;
    end
    
    % Separate into velocity/pressure
    n_dofs_u = u.nDofs;
    vars.u = fullx(1:n_dofs_u);
    vars.p = fullx(n_dofs_u+1:end);
    
    % Assign to FE functions
    n_node = round(length(vars.u) / u.ndimf);
    velfval = reshape(vars.u, u.ndimf, n_node)';
    u.setFValues(velfval)
    p.setFValues(vars.p)
    
    % Check convergence
    norm_dx = norm(dx) / (norm(x(free_dofs)) + eps);
    fprintf('Picard Iteration %d: Relative Norm dx = %.3e\n', iter, norm_dx);
    if norm_dx < tol
        fprintf('Converged after %d Picard iterations.\n', iter);
        break;
    end
end

% =========================================================================
% 7) EXTRACT AND PLOT THE RESULTS
% =========================================================================
uD  = dirichlet(:,3);
if ~isempty(dir_dofs)
    x(dir_dofs) = uD;
end
n_dofs_u = u.nDofs;
vars.u = x(1:n_dofs_u);
vars.p = x(n_dofs_u+1:end);

% Define final velocity/pressure
n_node = round(length(vars.u)/u.ndimf);
velfval = reshape(vars.u, u.ndimf, n_node)';
u.setFValues(velfval)
p.setFValues(vars.p)

u.plot();
p.plot();

% Plot velocity profiles at mid-sections
mid_x = 0.5;
mid_idx = find(abs(mesh.coord(:,1) - mid_x) < 1e-3);
figure;
plot(mesh.coord(mid_idx,2), velfval(mid_idx,1), 'LineWidth', 2);
xlabel('y-coordinate');
ylabel('u-velocity');
title('Velocity Profile at x = 0.5');
grid on;

mid_y = 0.5;
mid_idy = find(abs(mesh.coord(:,2) - mid_y) < 1e-3);
figure;
plot(mesh.coord(mid_idy,1), velfval(mid_idy,2), 'LineWidth', 2);
xlabel('x-coordinate');
ylabel('v-velocity');
title('Velocity Profile at y = 0.5');
grid on;
