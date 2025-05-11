clear
clc

% =========================================================================
% =========================================================================
%                      COMPRESSIBLE LINEAR ELASTICITY
% =========================================================================
% =========================================================================

% =========================================================================
% 1) CREATE MESH
% =========================================================================
mesh = UnitTriangleMesh(100,100);

% =========================================================================
% 2) CREATE FINITE ELEMENT FUNCTIONS (velocity u, pressure p)
% =========================================================================
u = LagrangianFunction.create(mesh, 2, 'P1');
n_dofs = u.nDofs;

% =========================================================================
% 3) SET UP BOUNDARY CONDITIONS (Dirichlet)
% =========================================================================
isLeft   = @(coor) (abs(coor(:,1) - min(coor(:,1))) < 1e-12);
isRight  = @(coor) (abs(coor(:,1) - max(coor(:,1))) < 1e-12);
isBottom = @(coor) (abs(coor(:,2) - min(coor(:,2))) < 1e-12);
isTop    = @(coor) (abs(coor(:,2) - max(coor(:,2))) < 1e-12);

dir.domain    = @(coor) isLeft(coor);
dir.direction = [1,2];
dir.value     = [0,0]; 
dirichlet = DirichletCondition(mesh, dir);

neu.domain = @(coor) isRight(coor);
neu.direction = [1,2];
neu.value     = [0,-1]; 
neumann = DirichletCondition(mesh, neu);

c.mesh = mesh;
c.dirichletFun = dirichlet;
c.pointloadFun = neumann;
c.periodicFun = [];
bc = BoundaryConditions(c);

E         = ConstantFunction.create(1,mesh);
nu        = ConstantFunction.create(1/3,mesh);
s.ptype   = 'ELASTIC';
s.pdim    = '2D';
s.nelem   = mesh.nelem;
s.mesh    = mesh;
s.young   = E;
s.poisson = nu;
s.type = 'ISOTROPIC';
s.ndim = mesh.ndim;
material = Material.create(s);

% =========================================================================
% 5) BUILD MATRICES
% =========================================================================
% 1) Diffusion (viscous) part
s1.type  = 'ElasticStiffnessMatrix';
s1.mesh  = mesh;
s1.test     = u;
s1.trial    = u;
s1.material = material;
s1.quadratureOrder = 2;
LHS1 = LHSintegrator.create(s1);
K = LHS1.compute();

% 2) Forces vector
s2.type     = 'Elastic';
s2.scale    = 'MACRO';
s2.dim.ndimf = mesh.ndim;
s2.dim.nnodes = mesh.nnodes;
s2.dim.ndofs = s2.dim.nnodes * s2.dim.ndimf;
s2.dim.nnodeElem = mesh.nnodeElem;
s2.dim.ndofsElem = s2.dim.ndimf * s2.dim.nnodeElem;
s2.BC       = bc;
s2.mesh     = mesh;
s2.material = material;
s2.globalConnec = mesh.connec;
RHSint = RHSintegrator.create(s2);
F = RHSint.compute();


% =========================================================================
% 6) SOLVE SYSTEM (Monolithic)
% =========================================================================

all_dofs  = 1:u.nDofs;
free = setdiff(all_dofs , bc.dirichlet_dofs);

U = zeros(u.nDofs,1);

maxIter = 2000000;
rtol = 3e-6;

K_r  = K(free,free);
F_r  = F(free);

tau = 0.32;

for k = 1:maxIter
    U_new = U;

    % -- gradient step on the FREE dofs only
    grad = K_r*U(free) - F_r;
    U_new(free) = U(free) - tau*grad;
    
    % -- proximal operator = impose Dirichlet values
    U_new(bc.dirichlet_dofs) = bc.dirichlet_vals;
        
    % -- convergence test
    res = norm(U_new-U) / max(1,norm(U));
    if res < rtol
        fprintf('Converged in %d iterations\n',k);
        break
    end

    if mod(k,100) == 0
        fprintf('Residual at iteration %i: %d\n',k, res);
    end

    U = U_new;
end

u.setFValues(reshape(U_new,2,[])');

% =========================================================================
% 7) PLOT RESULTS
% =========================================================================
u.plot();
