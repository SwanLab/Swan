% =========================================================================
% INCOMPRESSIBLE NAVIER–STOKES WITH PICARD + PROJECTION
% =========================================================================
clear; clc;

% CREATE MESH
mesh = UnitTriangleMesh(100,100);

% DEFINE FUNCTIONS
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

% BUILD CONSTANT MATRICES
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
dt = 1e-1;
nu = 1;
maxPicardIts = 20;
utol = 1e-7;
ptol = 1;
alpha = 0.7;

% initial
u_k = reshape(u.fValues', [],1);
p_k = reshape(p.fValues', [],1);
K_p(1,:)=0; K_p(:,1)=0; K_p(1,1)=1;

% TIME LOOP
for t = 1:time_steps
    % Picard initialize
    u_old = u_k;
    p_old = p_k;

    for pic = 1:maxPicardIts
        u.setFValues(reshape(u_old,2,[])');
        N = ConvForm.compute(u);

        A = (1/dt)*M - N + nu*K_u;
        LHS = [A, C; C', sparse(n_dir,n_dir)];
        RHS = [(1/dt)*M*u_k; dirichlet(:,3)];
        sol = LHS \ RHS;
        u_ast = sol(1:u_nd);

        rhs_p = D*u_ast/dt;
        rhs_p(1)=0;
        p_new = K_p \ rhs_p;

        x = [M C; C' sparse(n_dir,n_dir)] \ [(M*u_ast + dt * (G * p_new)); dirichlet(:,3)];
        u_corr = x(1:u.nDofs);
        u_new = u_old + alpha*(u_corr - u_old);

        velRes = norm(u_new-u_old)/(norm(u_new)+eps);
        divRes = norm(D*u_new);
        fprintf('t=%d pic=%d: uΔ=%.2e, div=%.2e\n',t,pic,velRes,divRes);
        if velRes<utol && divRes<ptol
            break;
        end

        u_old = u_new;
        p_old = p_new;
    end

    u_k = u_new;
    p_k = p_new;
    u.setFValues(reshape(u_k,2,[])');
    p.setFValues(reshape(p_k,1,[])');
    fprintf('Completed time step %d in %d Picard iterations.\n',t,pic);
end

u.plot();
p.plot();
