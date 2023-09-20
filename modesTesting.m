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
material = createMaterial(mesh,1);
% s.dim = '2D';
% bc = bcV;
nDimf=2;

dispFun=P1Function.create(mesh, nDimf);

dispFun=P1Function.create(mesh,1);


K    = computeStiffnessMatrix(mesh,material,dispFun);
Kred = bc.fullToReducedMatrix(K);

M= computeMassMatrix(mesh,dispFun);
Mred= bc.fullToReducedMatrix(M);
[basis,D]=eigs(M\K);


for i = 1:size(basis,2)
b = basis(:,i);
% b1 = bc.reducedToFullVector(b);
% bC{i} = reshape(b1,2,[])';
 bC{i} = b;

% sF.fValues = reshape(b1,2,[])';
% sF.mesh    = mesh;
% bF{i} =P1Function(sF);
% bF.plot
end

sM.mesh    = mesh;
sM.basis   = bC;
sM.fValues =[1 0 0 0 0 0];
modal = ModalFunction(sM);

p1FUNC = modal.project('P1');
% p1FUNC.plot
p1FUNC.print('prova')

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
file = 'CantileverBeam_Triangle_Linear_Fine';
file = 'Cantileverbeam_Quadrilateral_Bilinear';
a.fileName = file;
s = FemDataContainer(a);
mesh = s.mesh;


% % Generate coordinates
% x1 = linspace(0,2,20);
% x2 = linspace(0,1,20);
% % Create the grid
% [xv,yv] = meshgrid(x1,x2);
% % Triangulate the mesh to obtain coordinates and connectivities
% [F,basis] = mesh2tri(xv,yv,zeros(size(xv)),'x');
% 
% s.coord = basis(:,1:2);
% s.connec = F;
% mesh = Mesh(s);
end


function dim = getFunDims(mesh)
s.fValues = mesh.coord;
s.mesh = mesh;
disp = P1Function(s);
d.ndimf  = disp.ndimf;
d.ndimf  = 1;
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
s.type     = 'StiffnessMatrix';
s.mesh     = mesh;
% s.fun      = displacementFun;
s.test      = displacementFun;
s.trial      = displacementFun;
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