
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
isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-12);
isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1)))   < 1e-12);
isBottom = @(coor) (abs(coor(:,2) - min(coor(:,2)))   < 1e-12);
isTop    = @(coor) (abs(coor(:,2) - max(coor(:,2)))   < 1e-12);

dir_vel{1}.domain    = @(coor) isTop(coor);
dir_vel{1}.direction = [1,2];
dir_vel{1}.value     = [100,0];

dir_vel{2}.domain    = @(coor) (isBottom(coor) | isRight(coor) | isLeft(coor)) & not(isTop(coor));
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
sAF.fHandle = @(coor) [0.*coor,0.*coor];
sAF.ndimf   = 2;
sAF.mesh    = mesh;
forcesFormula = AnalyticalFunction(sAF);

% =========================================================================
% 5) BUILD MATRICES FOR THE LINEAR STOKES PROBLEM
% =========================================================================
% (a) Diffusion / stiffness matrix
s1.type  = 'StiffnessMatrix';
s1.mesh  = mesh;
s1.test  = u;
s1.trial = u;
LHS1 = LHSIntegrator.create(s1);
K = LHS1.compute();

% (b) Divergence matrix
s2.type = 'WeakDivergence';
s2.mesh = mesh;
s2.trial = p;
s2.test  = u;
LHS2 = LHSIntegrator.create(s2);
D = LHS2.compute();

% (c) Forces vector
s3.type          = 'Stokes';
s3.mesh          = mesh;
s3.velocityFun   = u;
s3.pressureFun   = p;
s3.forcesFormula = forcesFormula;
RHSint = RHSIntegrator.create(s3);
F = RHSint.integrate();

% (d) Nonlinear convective term
s4.type = 'NonLinearNS';
s4.mesh = mesh;
s4.trial = u;
s4.test = u;
LHS4 = LHSIntegrator.create(s4);

% C
n_dir = length(dir_dofs);
C = sparse(dir_dofs, 1:n_dir, 1, u.nDofs, n_dir);

% 0
O1 = zeros(p.nDofs);
O2 = zeros(p.nDofs, n_dir);
O3 = zeros(n_dir);
O4 = zeros(p.nDofs, 1);

RHS = [F; dirichlet(:,3)];
pFix = 1;
RHS(pFix + u.nDofs) = 0;

% =========================================================================
% 6) SOLVE SYSTEM
% =========================================================================
maxIter = 100;
tol = 1e-8;
x = zeros(length(RHS),1);
x_prev = zeros(length(RHS),1);
lambda = 1;

for iter = 1:maxIter

    N = LHS4.compute(u);
    
    LHS = [K+N D C; D' O1 O2; C' O2' O3];
    LHS(pFix + u.nDofs, :) = 0;
    LHS(:, pFix + u.nDofs) = 0;
    LHS(pFix + u.nDofs, pFix + u.nDofs) = 1;

    x = LHS\RHS;
    if ~isempty(dir_dofs)
        x(dir_dofs) = dirichlet(:,3);
    end
    
    vars.u = x(1:u.nDofs);
    vars.p = x(u.nDofs+1:end);
    
    % Assign to FE functions
    n_node = round(length(vars.u) / u.ndimf);
    velfval = reshape(vars.u, u.ndimf, n_node)';
    u.setFValues(velfval)
    p.setFValues(vars.p)
    
    % Check convergence
    norm_dx = norm(x - x_prev) / (norm(x) + eps);
    fprintf('Picard Iteration %d: Relative Norm dx = %.3e\n', iter, norm_dx);
    if norm_dx < tol
        fprintf('Converged after %d Picard iterations.\n', iter);
        break;
    end
    x_prev = x;
end

% =========================================================================
% 7) EXTRACT AND PLOT THE RESULTS
% =========================================================================
uD  = dirichlet(:,3);
if ~isempty(dir_dofs)
    x(dir_dofs) = uD;
end
vars.u = x(1:u.nDofs);
vars.p = x((1:p.nDofs) + u.nDofs);

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
