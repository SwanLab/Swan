function modesTesting


clc;clear all;close all

addpath(genpath(fileparts(mfilename('fullpath'))))

%file = 'test2d_triangle';
%a.fileName = file;
% s = FemDataContainer(a);


mesh = createMesh();

bcV = createRawBoundaryConditions(mesh);
bc = createBoundaryConditions(mesh,bcV);

% s.mesh = mesh;
% s.type = 'StiffnessMatrix';
% s.scale = 'MACRO';


quad = Quadrature.set(mesh.type);
quad.computeQuadrature('QUADRATIC');


material = createMaterial(mesh,quad.ngaus);
% s.dim = '2D';
% bc = bcV;
nDimf=2;

dispFun=P1Function.create(mesh, nDimf);


K    = computeStiffnessMatrix(mesh,material,dispFun);
Kred = bc.fullToReducedMatrix(K);

dim.ndimf  = dispFun.ndimf;
dim.nnodes = size(dispFun.fValues, 1);
dim.ndofs  = dim.nnodes*dim.ndimf;
dim.nnodeElem = mesh.nnodeElem; % should come from interp..
dim.ndofsElem = dim.nnodeElem*dim.ndimf;
c.dim=dim;
c.mesh=mesh;
c.BC = bc;
RHS    = RHSintegrator_ElasticMacro(c);
Fext = RHS.compute();
Fred = bc.fullToReducedVector(Fext);

M= computeMassMatrix(mesh,dispFun);
Mred= bc.fullToReducedMatrix(M);
% [basis,D]=eigs(M\K);
[basis,D]=eigs(Kred,8,'smallestabs');
psi = Kred*basis;

for i = 1:size(basis,2)
    b = basis(:,i);
    b1 = bc.reducedToFullVector(b);
    bC{i} = reshape(b1,2,[])';
    a = psi(:,i);
    a1=bc.reducedToFullVector(a);
    psiD{i}=reshape(a1,2,[])';
    %  bC{i} = b;

    % sF.fValues = reshape(b1,2,[])';
    % sF.mesh    = mesh;
    % bF{i} =P1Function(sF);
    % bF.plot
end
kbb=basis'*Kred*basis;
functionType = {'P1' , 'P1' , 'P1' , 'P1',  'P1',  'P1', 'P1', 'P1', 'P1', 'P1'};

sM.mesh    = mesh;
sM.basis   = bC;
sM.fValues =[1 1 1 1 1 1 1 1 1 1];
sM.functionType = {'P1' , 'P1' , 'P1' , 'P1',  'P1',  'P1', 'P1', 'P1', 'P1', 'P1'};
modal = ModalFunction(sM);
p1FUNC = modal.project('P1');
modal2=p1FUNC.project('ModalFunction',bC,functionType);


modalFun=ModalFunction.create(mesh,bC,functionType);

sL.material=material;
sL.test= modalFun;
sL.trial=modalFun;
sL.mesh=mesh;
sL.quadratureOrder = quad.order;
LHS=LHS_integratorStiffnessGlobal(sL);
lhs=LHS.compute();
LHSmass = LHS_integratorMassGlobal(sL);
lhsMass = LHSmass.compute();


%% create eifem matrices
bMesh = mesh.createBoundaryMesh();

rbDom = createRBfun(mesh);
rbBd = createRBfun(bMesh);
phib = createBoundaryDefFun(bMesh,bC,functionType);
psib = createBoundaryDefFun(bMesh,psiD,functionType);
nbound = size(bMesh,1);

for ibound=1:nbound
    sL.test  = phib{ibound};
    sL.trial = rbBd{ibound};
    sL.mesh  = bMesh{ibound}.mesh;
    sL.quadratureOrder = quad.order;
    lhs2{ibound} = LHS_integratorMassGlobal(sL);
    Hhat{ibound} = lhs2{ibound}.compute();
end

for ibound=1:nbound
    sL.test  = phib{ibound};
    sL.trial = psib{ibound};
    sL.mesh  = bMesh{ibound}.mesh;
    sL.quadratureOrder = quad.order;
    lhs2{ibound} = LHS_integratorMassGlobal(sL);
    H{ibound} = lhs2{ibound}.compute();
end

for ibound=1:nbound
    sL.test  = rbBd{ibound};
    sL.trial = rbBd{ibound};
    sL.mesh  = bMesh{ibound}.mesh;
    sL.quadratureOrder = quad.order;
    lhs2{ibound} = LHS_integratorMassGlobal(sL);
    G{ibound} = lhs2{ibound}.compute();
end

for ibound=1:nbound
    sL.test  = psib{ibound};
    sL.trial = rbBd{ibound};
    sL.mesh  = bMesh{ibound}.mesh;
    sL.quadratureOrder = quad.order;
    lhs2{ibound} = LHS_integratorMassGlobal(sL);
    Gtest{ibound} = lhs2{ibound}.compute(); %this should be 0 accoring joaquin's article
end

%%
% kprecon =

% refPoint = [0,0];
% RB=RigidBodyFunction.create(mesh,refPoint);
%
% sL.material=material;
% sL.test= RB;
% sL.trial=modalFun;
% sL.mesh=mesh;k
% sL.quadratureOrder = quad.order;
% LHS=LHS_integratorStiffnessGlobal(sL);
% lHS=LHS.compute();




% p1FUNC.plot
% p1FUNC.print('prova')

% modal.FEfun{1}.plot

%
% s.basis = bF;
% s.coef  = [2 3 -1];
% mF = ModalFunction(s); (%evaluate defined)
%
% m1 = mF.project('P1');
% m1.plot
%
% mF.plot()

[xCG,residualCG] = conjugateGradient(Kred,Fred);

%M = basis*inv(lhs)*basis';
%M = sparse(eye(size(xCG,1))./diag(Kred));

%M = basis*lhs*basis';

M2 = spdiags(1./sqrt(diag((Kred))),0,size(Kred,1),size(Kred,1));
L = ichol(Kred);
Apr = sparse(diag(diag(Kred))\Kred);
bpr = sparse(diag(diag(Kred))\Fred);
M=basis*(lhs\basis');
Mb{1}=basis';
Mb{2}= inv(lhs);
Mb{3}=basis;
% Mb{1}=eye(size(Kred,1));
% M=M+M2;
%[xCG,residualCG] = conjugateGradient(Apr,bpr);
% xCG2 = pcg(Apr,bpr,1E-6,1000);
% [xCG2,residualCG2] = conjugateGradient(Kred,Fred,M2);

[xCG2,residualCG2] = preconditionedConjugateGradient(Kred,Fred,Mb);

x= Kred\Fred;
[X1,FLAG1,RELRES1,iter1] = pcg(Kred,Fred,1E-6,1000,L,L');
[X2,FLAG2,RELRES2,iter2] = pcg(Kred,Fred,1E-6,1000,M2,M2');
[X3,FLAG3,RELRES3,iter3] = pcg(Kred,Fred,1E-6,1000);
s.mesh = mesh;
s.type = 'ELASTIC';
s.scale = 'MACRO';
s.material = createMaterial(mesh,1);
s.dim = '2D';
s.bc = bcV;
s.nDimf=2;


fem = FEM.create(s);
fem.solve();

fem.print('DD')



% figure(1)
% fem.uFun.plot()
% figure(2)
% fem.stressFun.plot()
% figure(3)
% fem.strainFun.plot()


% FEM
%solution plot

% eigmodes LHS

% ModalFunction build


% Plot an example with given values


end

function mesh = createMesh()
% file = 'CantileverBeam_Triangle_Linear';
% % file = 'Cantileverbeam_Quadrilateral_Bilinear';
% a.fileName = file;
% s = FemDataContainer(a);
% mesh = s.mesh;


% Generate coordinates
x1 = linspace(0,2,20);
x2 = linspace(0,1,20);
% Create the grid
[xv,yv] = meshgrid(x1,x2);
% Triangulate the mesh to obtain coordinates and connectivities
[F,basis] = mesh2tri(xv,yv,zeros(size(xv)),'x');

s.coord = basis(:,1:2);
s.connec = F;
mesh = Mesh(s);
end


function dim = getFunDims(mesh)
s.fValues = mesh.coord;
s.mesh = mesh;
disp = P1Function(s);
d.ndimf  = disp.ndimf;
d.nnodes = size(disp.fValues, 1);
d.ndofs  = d.nnodes*d.ndimf;
d.nnodeElem = mesh.nnodeElem; % should come from interp..
d.ndofsElem = d.nnodeElem*d.ndimf;
dim = d;
end

function bc = createBoundaryConditions(mesh,bcV)
dim = getFunDims(mesh);
bcV.ndimf = dim.ndimf;
bcV.ndofs = dim.ndofs;
s.mesh  = mesh;
s.scale = 'MACRO';
s.bc    = {bcV};
s.ndofs = dim.ndofs;
bc = BoundaryConditions(s);
bc.compute();
end



function bc = createRawBoundaryConditions(mesh)
dirichletNodes = abs(mesh.coord(:,1)-0) < 1e-12;
rightSide  = max(mesh.coord(:,1));
isInRight = abs(mesh.coord(:,1)-rightSide)< 1e-12;
isInMiddleEdge = abs(mesh.coord(:,2)-0.5) < 0.1;
forceNodes = isInRight & isInMiddleEdge;
forceNodes = isInRight;
nodes = 1:mesh.nnodes;
bcDir = [nodes(dirichletNodes)';nodes(dirichletNodes)'];
% bcDir = [nodes(dirichletNodes)'];

nodesdir=size(nodes(dirichletNodes),2);
bcDir(1:nodesdir,end+1) = 1;
bcDir(nodesdir+1:end,end) = 2;
bcDir(:,end+1)=0;
bc.dirichlet = bcDir;
bc.pointload(:,1) = nodes(forceNodes);
bc.pointload(:,2) = 2;
bc.pointload(:,3) = -1;
end


function material = createMaterial(mesh,ngaus)
I = ones(mesh.nelem,ngaus);
s.ptype = 'ELASTIC';
s.pdim  = '2D';
s.nelem = mesh.nelem;
s.mesh  = mesh;
s.kappa = .9107*I;
s.mu    = .3446*I;
mat = Material.create(s);
mat.compute(s);
material = mat;
end


function k= computeStiffnessMatrix(mesh,material,displacementFun)
s.type     = 'ElasticStiffnessMatrix';
s.mesh     = mesh;
s.fun      = displacementFun;
% s.test      = displacementFun;
% s.trial      = displacementFun;
s.material = material;
lhs = LHSintegrator.create(s);
k=lhs.compute();
end

function m=computeMassMatrix(mesh,displacementFun)
s.type     = 'MassMatrix';
s.mesh     = mesh;
s.test      = displacementFun;
s.trial      = displacementFun;
s.quadratureOrder = 'QUADRATIC';
lhs = LHSintegrator.create(s);
m=lhs.compute();
end

function [x,residual] = conjugateGradient(A,B,M)
if nargin == 3
    LHS=(M*A*M');
    RHS=M*B;
else
    LHS=A;
    RHS=B;
    M = speye(size(B,1));
end
tol = 1e-6;
n = length(RHS);
x = zeros(n,1);
r = RHS - LHS * x;
p = r;
rsold = r' * r;
iter = 0;

hasNotConverged = true;

while hasNotConverged
    Ap = LHS * p;
    alpha = rsold / (p' * Ap);
    x = x + alpha * p;
    r = r - alpha * Ap;
    rsnew = r' * r;

    %hasNotConverged = sqrt(rsnew) > tol;
    hasNotConverged = norm(LHS*x - RHS) > tol;

    p = r + (rsnew / rsold) * p;
    rsold = rsnew;
    iter = iter + 1;
    residual(iter) = norm(LHS*x - RHS); %Ax - b
end
x=(M)*x;
end


function [x,residual] = preconditionedConjugateGradient(A,B,M)
tol = 1e-6;
n = length(B);
x = zeros(n,1);
r = B - A * x;
z = matVecProd(M,r);
z=r-z;
p = z;
rzold = r' * z;
iter = 0;

hasNotConverged = true;

while hasNotConverged
    Ap = A * p;
    alpha = rzold / (p' * Ap);
    x = x + alpha * p;
    r = r - alpha * Ap;
    z = matVecProd(M,r);
    z = r-z;
    rznew = r' * z;

    %hasNotConverged = sqrt(rsnew) > tol;
    hasNotConverged = norm(r) > tol;

    p = z + (rznew / rzold) * p;
    rzold = rznew;
    iter = iter + 1;
    residual(iter) = norm(r); %Ax - b
end
% x=(M)*x;
end

function v = matVecProd(M,x)
    for i = 1:size(M,2)
        x=M{i}*x;
    end
    v = x;
end


function rbFun = createRBfun(mesh)
%domain
if size(mesh,1)==1
    centroid = [sum(mesh.coord(:,1))/mesh.nnodes,sum(mesh.coord(:,2))/mesh.nnodes];
    rbFun = RigidBodyFunction.create(mesh,centroid);
else
    %boundary
    for i=1:size(mesh,1)
        meshb = mesh{i}.mesh;
        centroid = [sum(meshb.coord(:,1))/meshb.nnodes,sum(meshb.coord(:,2))/meshb.nnodes];
        rbFun{i} = RigidBodyFunction.create(meshb,centroid);
    end
end
end


function defFun = createBoundaryDefFun(bmesh,basis,functionType)
for imesh=1:size(bmesh,1)
    nodes = unique(bmesh{imesh}.globalConnec);
    %         nnode= length(nodes);
    %         for inode=1:nnode
    %             dof = 2*(nodes(inode)-1)+1;
    %             dofs(dof)=dof;
    %             dofs(dof+1)=dof+1;
    %         end
    %         dofs=dofs(dofs>0);
    for ibasis=1:length(basis)
        basisb{ibasis} = basis{ibasis}(nodes,:);
    end

    meshb = bmesh{imesh}.mesh;
    defFun{imesh} = ModalFunction.create(meshb,basisb,functionType);
    clear dofs;
end
end



%
% function a = ichol(a)
% n = size(a,1);
%
% for k = 1:n
%     a(k,k) = sqrt(a(k,k));
%     for i = (k+1):n
%         if (a(i,k) ~= 0)
%             a(i,k) = a(i,k)/a(k,k);
%             end
%             end
%             for j = (k+1):n
%                 for i = j:n
%                     if (a(i,j) ~= 0)
%                         a(i,j) = a(i,j) - a(i,k)*a(j,k);
%                     end
%                 end
%             end
%         end
%
%         for i = 1:n
%             for j = i+1:n
%                 a(i,j) = 0;
%             end
%         end
%     end
% end
% end