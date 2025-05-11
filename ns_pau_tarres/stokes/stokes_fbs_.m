clear
clc

% =========================================================================
% 1) CREATE MESH
% =========================================================================
mesh = UnitTriangleMesh(70,70);

% =========================================================================
% 2) CREATE FINITE ELEMENT FUNCTIONS (velocity u, pressure p)
% =========================================================================
u = LagrangianFunction.create(mesh, 2, 'P2');  % velocity
p = LagrangianFunction.create(mesh, 1, 'P1');  % pressure
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
dir_vel{1}.value     = [1,0];    % velocity on the top boundary

dir_vel{2}.domain    = @(coor) isBottom(coor) | isRight(coor) | isLeft(coor);
dir_vel{2}.direction = [1,2];
dir_vel{2}.value     = [0,0];    % velocity = 0 elsewhere

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
% (a) Diffusion / stiffness matrix
s1.type  = 'StiffnessMatrix';
s1.mesh  = mesh;
s1.test  = u;
s1.trial = u;
LHS1 = LHSintegrator.create(s1);
K = LHS1.compute();

% (b) Divergence matrix
s2.type  = 'WeakDivergence';
s2.mesh  = mesh;
s2.trial = p;
s2.test  = u;
LHS2 = LHSintegrator.create(s2);
D = LHS2.compute();

% (c) Forces vector
s3.type          = 'Stokes';
s3.mesh          = mesh;
s3.velocityFun   = u;
s3.pressureFun   = p;
s3.forcesFormula = forcesFormula;
RHSint = RHSintegrator.create(s3);
f = RHSint.integrate();
f_u = f(1:n_dofs_u);
f_p = f((n_dofs_u+1):n_dofs);


% =========================================================================
% 6) FORWARD-BACKWARD ITERATION
% =========================================================================
maxIter = 1000;
tol     = 1e-6;

% Choose a step size λ < 2/L. We approximate the largest eigenvalue of LHS:
lambda = 0.005;

u_k1 = zeros(n_dofs_u, 1);
p_k1 = zeros(n_dofs_p, 1);

for iter = 1:maxIter
    u_k = u_k1;
    p_k = p_k1;

    % (1) FORWARD STEP
    grad_F = K*u_k - f_u;
    u_half = u_k - lambda*grad_F;
    
    % (2) BACKWARD STEP (prox)
    p_k1 = (D'*D)\( D'*u_half / lambda );
    u_k1 = u_half - lambda*( D * p_k1 );

    u_k1(dir_dofs) = dirichlet(:,3);


    step_norm = norm(u_k1 - u_k) / (norm(u_k) + eps);
    if step_norm < tol
        fprintf('FBS converged at iteration %d, rel. step = %.3e\n',iter,step_norm);
        break;
    end

    if mod(iter,50)==0
        fprintf('Iter: %i\t StepNorm: %d\n', iter,step_norm)
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
