% =========================================================================
% INCOMPRESSIBLE NAVIER–STOKES WITH PICARD + PROJECTION
% =========================================================================
clear; clc;

% CREATE MESH
mesh = UnitTriangleMesh(100,100);

% DEFINE FE SPACES
u     = LagrangianFunction.create(mesh, 2, 'P2');
p     = LagrangianFunction.create(mesh, 1, 'P1');
u_nd  = u.nDofs;
p_nd  = p.nDofs;

% BOUNDARY CONDITIONS for velocity (Dirichlet)
isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1)))   < 1e-12);
isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1)))   < 1e-12);
isBottom = @(coor) (abs(coor(:,2) - min(coor(:,2)))   < 1e-12);
isTop    = @(coor) (abs(coor(:,2) - max(coor(:,2)))   < 1e-12);

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
    nodes = 1 + (dirDofs(2:2:end) - 2) / u.ndimf;
    iNod = repelem(nodes, 2)';
    dofType = repmat((1:2)', length(nodes), 1);
    valmat = repmat(dir_vel{i}.value', length(nodes), 1);
    dirichlet = [dirichlet; [iNod', dofType, valmat]];
    dir_dofs = [dir_dofs; dirDofs(:)];
end

% FORCES
sf.mesh    = mesh;
sf.ndimf   = 2;
sf.fHandle = @(coor) [0*coor,0*coor];
forces = AnalyticalFunction(sf);

% BUILD MATRICES
group.type  = 'MassMatrix';
group.mesh  = mesh;
group.test = u;
group.trial = u;
M = LHSIntegrator.create(group).compute();

group.type = 'StiffnessMatrix';
group.test = u;
group.trial = u;
K_u = LHSIntegrator.create(group).compute();

group.type = 'StiffnessMatrix';
group.test = p;
group.trial = p;
K_p = LHSIntegrator.create(group).compute();

group.type = 'NonLinearNS';
group.trial = u;
group.test = u;
ConvForm = LHSIntegrator.create(group);

group.type = 'WeakDivergence';
group.mesh = mesh;
group.test = u;
group.trial = p;
D_weak = LHSIntegrator.create(group).compute();
D = D_weak'; % div
G = -D'; % gradient

n_dir = length(dir_dofs);
C = sparse(dir_dofs, 1:n_dir, 1, u.nDofs, n_dir);

% TIME-STEPPING SETUP
time_steps = 100;
dt = 1e-4;
nu = 1;
maxPicardIts  = 300;
utol = 1e-10;
alpha = 0.7;

u_k = reshape(u.fValues', [],1);
K_p(1,:)=0; K_p(:,1)=0; K_p(1,1)=1;

for pic = 1:maxPicardIts
    u.setFValues(reshape(u_k,2,[])');
    N = ConvForm.compute(u);

    A = (1/dt)*M - N + nu*K_u;
    LHS = [A, C; C', sparse(n_dir,n_dir)];
    RHS = [(1/dt)*M*u_k; dirichlet(:,3)];
    sol = LHS \ RHS;
    u_ast = sol(1:u_nd);

    rhs_p = D*u_ast/dt;
    rhs_p(1)=0;
    p_k1 = K_p \ rhs_p;

    x = [M C; C' sparse(n_dir,n_dir)] \ [(M*u_ast + dt * (G * p_k1)); dirichlet(:,3)];
    u_k1 = x(1:u.nDofs);
    u_k1 = u_k + alpha*(u_k1 - u_k);

    velRes = norm(u_k1-u_k)/(norm(u_k1)+eps);
    divRes = norm(D*u_k1);
    fprintf('it=%d: uΔ=%.2e, div=%.2e\n',pic,velRes,divRes);
    if velRes<utol
        break;
    end

    u_k = u_k1;
end

u_k = u_k1;
p_k = p_k1;
u.setFValues(reshape(u_k,2,[])');
p.setFValues(reshape(p_k,1,[])');

% PLOT RESULTS
u.plot();
p.plot();
