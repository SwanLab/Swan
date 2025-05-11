clear
clc

% =========================================================================
% 1) CREATE MESH
% =========================================================================
mesh = UnitTriangleMesh(30,30);

% =========================================================================
% 2) CREATE FINITE ELEMENT FUNCTIONS (velocity u, pressure p)
% =========================================================================
u = LagrangianFunction.create(mesh, 2, 'P2');
p = LagrangianFunction.create(mesh, 1, 'P1');
n_dofs_u = u.nDofs;
n_dofs_p = p.nDofs;
n_dofs   = n_dofs_u + n_dofs_p;

% =========================================================================
% 3) SET UP BOUNDARY CONDITIONS (Dirichlet)
% =========================================================================
isLeft   = @(coor) abs(coor(:,1) - min(coor(:,1)))   < 1e-12;
isRight  = @(coor) abs(coor(:,1) - max(coor(:,1)))   < 1e-12;
isBottom = @(coor) abs(coor(:,2) - min(coor(:,2)))   < 1e-12;
isTop    = @(coor) abs(coor(:,2) - max(coor(:,2)))   < 1e-12;

dir_vel{1}.domain    = @(coor) isTop(coor);
dir_vel{1}.direction = [1,2];
dir_vel{1}.value     = [1,0];

dir_vel{2}.domain    = @(coor) (isBottom(coor) | isRight(coor) | isLeft(coor)) & not(isTop(coor));
dir_vel{2}.direction = [1,2];
dir_vel{2}.value     = [0,0]; 

dirichlet = [];
dir_dofs = [];
for i = 1:length(dir_vel)
    dirDofs = u.getDofsFromCondition(dir_vel{i}.domain);
    dofs   = 1 + (dirDofs(2:2:end) - 2) / u.ndimf;
    iNod    = repelem(dofs, 2)';
    dofType = repmat((1:2)', length(dofs), 1);
    valmat  = repmat(dir_vel{i}.value', length(dofs), 1);
    dirichlet  = [dirichlet; [iNod', dofType, valmat]];
    dir_dofs   = [dir_dofs; dirDofs(:)];
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
s1.type  = 'MassMatrix';
s1.mesh  = mesh;
s1.test  = u;
s1.trial = u;
LHS1 = LHSIntegrator.create(s1);
M = LHS1.compute();

s2.type  = 'StiffnessMatrix';
s2.mesh  = mesh;
s2.test  = u;
s2.trial = u;
LHS2 = LHSIntegrator.create(s2);
K = LHS2.compute();

s3.type  = 'WeakDivergence';
s3.mesh  = mesh;
s3.trial = p;
s3.test  = u;
LHS3 = LHSIntegrator.create(s3);
B = LHS3.compute();

s4.type = 'Stokes';
s4.mesh = mesh;
s4.velocityFun = u;
s4.pressureFun = p;
s4.forcesFormula = forcesFormula;
RHSint = RHSIntegrator.create(s4);
F = RHSint.integrate();
F = F(1:n_dofs_u);

O = zeros(n_dofs_p,n_dofs_p);
o = zeros(n_dofs_p, 1);

LHS = [M B; B' O];

pFix = 1;
LHS(pFix + u.nDofs, :) = 0;
LHS(:, pFix + u.nDofs) = 0;
LHS(pFix + u.nDofs, pFix + u.nDofs) = 1;
RHS(pFix + u.nDofs) = 0;


% =========================================================================
% 6) FORWARD-BACKWARD ITERATION
% =========================================================================
maxIter = 200000;
tol     = 1e-6;

lambda = 0.15;

u_k1 = zeros(n_dofs_u, 1);
p_k1 = zeros(n_dofs_p, 1);

for iter = 1:maxIter
    u_k = u_k1;
    p_k = p_k1;
    u_k(dir_dofs) = dirichlet(:,3);
    u_ast = u_k - lambda*(K*u_k + F);
    rhs = [M*u_ast; o];
    vars = LHS \ rhs;
    u_k1 = vars(1:n_dofs_u);
    p_k1 = vars(1+n_dofs_u:end);

    % ---------- convergence monitors ----------
    step = norm(u_k1-u_k) / (norm(u_k) + eps);

    % ---------- stop tests ----------
    if step < 1e-6
        fprintf('Converged in %d its: step %.1e\n', iter, step);
        break
    end

    if mod(iter,10)==0
        fprintf('it %4d | step %.3e \n', iter, step);
    end
end

% =========================================================================
% 7) EXTRACT AND PLOT THE RESULTS
% =========================================================================
% Store velocity values in a matrix for post-processing
nu    = u.ndimf;
nNode = round(length(u_k1)/nu);
nodes = 1:nNode;
velfval = zeros(nNode, nu);
for idim = 1:nu
    dofs = nu*(nodes-1) + idim;
    velfval(:,idim) = u_k1(dofs);
end
u.setFValues(velfval);
p.setFValues(p_k1);

% Plot velocity and pressure solutions
u.plot();
p.plot();

% Example velocity profile at x=0.5
mid_x = 0.5;
mid_idx = find(abs(mesh.coord(:,1) - mid_x) < 1e-3);
figure;
plot(mesh.coord(mid_idx,2), velfval(mid_idx,1), 'LineWidth', 2);
xlabel('y-coordinate');
ylabel('u-velocity');
title('Velocity Profile at x = 0.5');
grid on;

% Example velocity profile at y=0.5
mid_y = 0.5;
mid_idy = find(abs(mesh.coord(:,2) - mid_y) < 1e-3);
figure;
plot(mesh.coord(mid_idy,1), velfval(mid_idy,2), 'LineWidth', 2);
xlabel('x-coordinate');
ylabel('v-velocity');
title('Velocity Profile at y = 0.5');
grid on;
