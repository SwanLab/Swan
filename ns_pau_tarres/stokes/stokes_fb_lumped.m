%% 0) HOUSEKEEPING
clear; clc;

% =========================================================================
% 1) CREATE MESH
% =========================================================================
mesh = UnitTriangleMesh(100,100);

% =========================================================================
% 2) CREATE FINITE-ELEMENT FUNCTIONS (velocity u, pressure p)
% =========================================================================
u = LagrangianFunction.create(mesh, 2, 'P2');
p = LagrangianFunction.create(mesh, 1, 'P1');
n_dofs_u = u.nDofs;
n_dofs_p = p.nDofs;
n_dofs   = n_dofs_u + n_dofs_p;

% =========================================================================
% 3) DIRICHLET BOUNDARY CONDITIONS
% =========================================================================
isLeft   = @(coor) abs(coor(:,1) - min(coor(:,1)))   < 1e-12;
isRight  = @(coor) abs(coor(:,1) - max(coor(:,1)))   < 1e-12;
isBottom = @(coor) abs(coor(:,2) - min(coor(:,2)))   < 1e-12;
isTop    = @(coor) abs(coor(:,2) - max(coor(:,2)))   < 1e-12;

dir_vel{1}.domain    = @(coor) isTop(coor);
dir_vel{1}.direction = [1,2];
dir_vel{1}.value     = [1,0];

dir_vel{2}.domain    = @(coor) (isBottom(coor) | isRight(coor) | isLeft(coor)) & ~isTop(coor);
dir_vel{2}.direction = [1,2];
dir_vel{2}.value     = [0,0]; 

dirichlet = [];
dir_dofs  = [];
for i = 1:length(dir_vel)
    dirDofs = u.getDofsFromCondition(dir_vel{i}.domain);
    dofs    = 1 + (dirDofs(2:2:end) - 2) / u.ndimf;
    iNod    = repelem(dofs, 2)';
    dofType = repmat((1:2)', length(dofs), 1);
    valmat  = repmat(dir_vel{i}.value', length(dofs), 1);
    dirichlet = [dirichlet; [iNod', dofType, valmat]];
    dir_dofs  = [dir_dofs; dirDofs(:)];
end

% =========================================================================
% 4) APPLIED BODY FORCES
% =========================================================================
sAF.fHandle = @(coor) [0.*coor, 0.*coor];
sAF.ndimf   = 2;
sAF.mesh    = mesh;
forcesFormula = AnalyticalFunction(sAF);

% =========================================================================
% 5) BUILDING MATRICES FOR THE LINEAR STOKES PROBLEM
% =========================================================================
% -- Mass matrix (velocity)
s1.type = 'MassMatrix'; s1.mesh = mesh; s1.test = u; s1.trial = u;
M  = LHSIntegrator.create(s1).compute();
M_lumped = spdiags(sum(M,2), 0, size(M,1), size(M,2));   % diagonal

% -- Stiffness matrix (velocity)
s2.type = 'StiffnessMatrix'; s2.mesh = mesh; s2.test = u; s2.trial = u;
K = LHSIntegrator.create(s2).compute();

% -- Divergence operator  B  (size n_dofs_u × n_dofs_p)
s3.type = 'WeakDivergence'; s3.mesh = mesh; s3.trial = p; s3.test = u;
B = LHSIntegrator.create(s3).compute();

% -- Load vector  F  (velocity part only)
s4.type = 'Stokes'; s4.mesh = mesh; s4.velocityFun = u;
s4.pressureFun = p; s4.forcesFormula = forcesFormula;
F = RHSIntegrator.create(s4).integrate();
F = F(1:n_dofs_u);

% -------------------------------------------------------------------------
% PRE-COMPUTE bits for the matrix-free Schur complement
% -------------------------------------------------------------------------
pFix = 1;                               % pressure anchor DOF
M_L_inv = 1 ./ full(diag(M_lumped));    % vector (n_dofs_u × 1)

% handle that applies  S = B' * M_L_inv * B   without assembling S
applyS = @(pvec) localSchur(pvec,B,M_L_inv,pFix);

% Jacobi preconditioner   P ≈ diag(S)
diagS = sum( (B.^2) .* M_L_inv , 1 )';   % column vector (n_dofs_p × 1)
diagS(pFix)  = 1;                               % keep pin row SPD
Minv = @(x) x ./ diagS;                 % function-handle M⁻¹x

% =========================================================================
% 6) FORWARD–BACKWARD (PROJECTION) ITERATION
% =========================================================================
maxIter  = 2e5;   tol = 1e-6;   lambda = 0.00001;
u_k1     = zeros(n_dofs_u,1);
p_k1     = zeros(n_dofs_p,1);

for iter = 1:maxIter
    % ---------------- velocity predictor ----------------
    u_k = u_k1;
    u_k(dir_dofs) = dirichlet(:,3);
    u_ast = u_k - lambda*(K*u_k - F);

    % ---------------- pressure solve  (Schur system) ----------------
    rhs_p       = B' * u_ast;
    rhs_p(pFix) = 0;                              % honour pinning

    tol_p   = 1e-8;
    max_pcg = 200;
    [p_k1,flag,rel,its] = pcg(applyS, rhs_p, tol_p, max_pcg, Minv, [], p_k1);

    if flag~=0
        warning('PCG did not converge: flag %d, relres %.2e, it %d',flag,rel,its);
    end

    % ---------------- velocity corrector ----------------
    u_k1 = u_ast - M_L_inv .* (B * p_k1);        % diag solve
    u_k1(dir_dofs) = dirichlet(:,3);             % re-enforce BC

    % ---------------- convergence check ----------------
    step = norm(u_k1-u_k) / (norm(u_k)+eps);
    if step < tol
        fprintf('Converged in %d iterations, step = %.2e\n',iter,step);
        break
    end
    if mod(iter,10)==0
        fprintf('it %5d | step %.3e | PCG its %3d\n',iter,step,its);
    end
end

% =========================================================================
% 7) POST-PROCESSING & PLOTS
% =========================================================================
nu     = u.ndimf;
nNode  = length(u_k1)/nu;
nodes  = 1:nNode;
velfval = zeros(nNode, nu);
for idim = 1:nu
    dofs = nu*(nodes-1) + idim;
    velfval(:,idim) = u_k1(dofs);
end
u.setFValues(velfval);
p.setFValues(p_k1);

u.plot();
p.plot();

% Example velocity profile at x = 0.5
mid_x   = 0.5;
mid_idx = find(abs(mesh.coord(:,1) - mid_x) < 1e-3);
figure;
plot(mesh.coord(mid_idx,2), velfval(mid_idx,1), 'LineWidth', 2);
xlabel('y'); ylabel('u-velocity'); title('u(y) at x = 0.5'); grid on;

% Example velocity profile at y = 0.5
mid_y   = 0.5;
mid_idy = find(abs(mesh.coord(:,2) - mid_y) < 1e-3);
figure;
plot(mesh.coord(mid_idy,1), velfval(mid_idy,2), 'LineWidth', 2);
xlabel('x'); ylabel('v-velocity'); title('v(x) at y = 0.5'); grid on;

% -------------------------------------------------------------------------
%  LOCAL FUNCTION: matrix-free Schur-complement operator
% -------------------------------------------------------------------------
function y = localSchur(pvec,B,M_L_inv,pFix)
    y        = B' * (M_L_inv .* (B * pvec));  %  Bᵀ M_L⁻¹ B p
    y(pFix)  = pvec(pFix);                    % enforce pin row (SPD)
end
