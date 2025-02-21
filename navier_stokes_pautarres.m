clear
clc

% CREATE MESH
mesh = UnitQuadMesh(1000,1000);

% CREATE FUNCTIONS
u = LagrangianFunction.create(mesh, 2, 'P2');
p = LagrangianFunction.create(mesh, 1, 'P1');
n_dofs = u.nDofs + p.nDofs;

% CREATE BOUNDARY CONDITIONS
isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-12);
isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1)))   < 1e-12);
isBottom = @(coor) (abs(coor(:,2) - min(coor(:,2)))   < 1e-12);
isTop    = @(coor) (abs(coor(:,2) - max(coor(:,2)))   < 1e-12);

dir_vel{1}.domain    = @(coor) isTop(coor);
dir_vel{1}.direction = [1,2];
dir_vel{1}.value     = [1,0];

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

% DEFINE APPLIED FORCES
sAF.fHandle = @(coor) [0.*coor,0.*coor];
sAF.ndimf   = 2;
sAF.mesh    = mesh;
forcesFormula = AnalyticalFunction(sAF);

% CREATE LEFT-HAND SIDE
s1.type  = 'StiffnessMatrix';
s1.mesh  = mesh;
s1.test  = u;
s1.trial = u;
LHS1 = LHSintegrator.create(s1);
K = LHS1.compute();

s2.type = 'WeakDivergence';
s2.mesh = mesh;
s2.trial = p;
s2.test  = u;
LHS2 = LHSintegrator.create(s2);
D = LHS2.compute();

% CREATE NONLINEAR CONVECTIVE TERM
s3.type = 'NonLinearNS';
s3.mesh = mesh;
s3.trial = u;
s3.test = u;
s3.velocityField = u;
LHS3 = LHSintegrator.create(s3);
C = LHS3.compute();


sz = size(D, 2);
BB = sparse(sz,sz);

LHS = [K+C D; D', BB];

% CREATE RIGHT-HAND SIDE
s5.type          = 'Stokes';
s5.mesh          = mesh;
s5.velocityFun   = u;
s5.pressureFun   = p;
s5.forcesFormula = forcesFormula;
RHSint = RHSintegrator.create(s5);
F = RHSint.integrate();
uD = dirichlet(:,3);
R  = -LHS(:,dir_dofs)*uD;
RHS = F + R;


% SOLVE PROBLEM
free_dofs_plus = setdiff(1:n_dofs,dir_dofs);
LHSr = LHS(free_dofs_plus,free_dofs_plus);
RHSr = RHS(free_dofs_plus);
x = LHSr\RHSr;





% SOLVE PROBLEM USING NEWTON-RAPHSON
maxIter = 30;
tol = 1e-6;
x = zeros(n_dofs,1);
free_dofs = setdiff(1:n_dofs, dir_dofs);

for iter = 1:maxIter
    u_vals = x(1:u.nDofs);
    u.setFValues(reshape(u_vals, [], u.ndimf));

    s3.velocityField = u;
    LHS3 = LHSintegrator.create(s3);
    C = LHS3.compute();
    LHS = [K + C, D; D', BB];

    R  = -LHS(:,dir_dofs) * uD;
    RHS = F + R - LHS * x;

    LHSr = LHS(free_dofs, free_dofs);
    RHSr = RHS(free_dofs);
    dx = LHSr\RHSr;

    x(free_dofs) = x(free_dofs) + dx;

    norm_dx = norm(dx) / (norm(x(free_dofs)) + eps);
    fprintf('Iteration %d: Relative Norm dx = %.3e\n', iter, norm_dx);
    if norm_dx < tol
        fprintf('Converged after %d iterations.\n', iter);
        break;
    end
end





% ADD DDIRICHLET BOUNDARY CONDITIONS
uD  = dirichlet(:,3);
nsteps = length(x(1,:));
uD = repmat(uD,1,nsteps);
free_dofs = setdiff(1:(n_dofs),dir_dofs);
fullx = x;
if ~isempty(dir_dofs)
    fullx(dir_dofs,:) = uD;
end

% SEPARATE VARIABLES FVALUES
n_dofs_u = u.nDofs;
vars.u = fullx(1:n_dofs_u,:);
vars.p = fullx(n_dofs_u+1:end,:);

% DEFINE VARIABLES
nu = u.ndimf;
nnode = round(length(vars.u)/nu);
nodes = 1:nnode;
velfval = zeros(nnode,nu);
for idim = 1:nu
    dofs = nu*(nodes-1)+idim;
    velfval(:,idim) = vars.u(dofs, end);
end
u.setFValues(velfval)
p.setFValues(vars.p(:,end))

% PLOT RESULTS
u.plot();
p.plot();

