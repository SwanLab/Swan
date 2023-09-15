function modesTesting
% file = 'test2d_triangle';
% a.fileName = file;
% s = FemDataContainer(a);
% mesh = s.mesh;

clc;clear;close all

addpath(genpath(fileparts(mfilename('fullpath'))))

%file = 'test2d_triangle';
%a.fileName = file;
% s = FemDataContainer(a);



% Generate coordinates
x1 = linspace(0,2,50);
x2 = linspace(0,1,25);
% Create the grid
[xv,yv] = meshgrid(x1,x2);
% Triangulate the mesh to obtain coordinates and connectivities
[F,basis] = mesh2tri(xv,yv,zeros(size(xv)),'x');

s.coord = basis(:,1:2);
s.connec = F;
mesh = Mesh(s);
bcV = createRawBoundaryConditions(mesh);
bc = createBoundaryConditions(mesh,bcV);

s.mesh = mesh;
s.type = 'ELASTIC';
s.scale = 'MACRO';
s.material = createMaterial(mesh,1);
s.dim = '2D';
s.bc = bcV;
s.nDimf=2;

dispFun=P1Function.create(s.mesh, s.nDimf);

K    = computeStiffnessMatrix(s.mesh,s.material,dispFun);
Kred = bc.fullToReducedMatrix(K);
[basis,D]=eigs((Kred));


for i = 1:3
b = basis(:,i);
b1 = bc.reducedToFullVector(b);

sF.fValues = reshape(b1,2,[])';
sF.mesh    = mesh;
bF{i} =P1Function(sF);
bF.plot
end


s.basis = bF;
s.coef  = [2 3 -1];
mF = ModalFunction(s); (%evaluate defined)

m1 = mF.project('P1');
m1.plot

mF.plot()


fem = FEM.create(s);
fem.solve();

fem.print('DD','Paraview')

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
nodes = 1:mesh.nnodes;
bc.dirichlet = nodes(dirichletNodes);
bc.pointload(:,1) = nodes(forceNodes);
bc.pointload(:,2) = 2;
bc.pointload(:,3) = -1;
end

function const = constDof(dirichlet,ndimf)
nnode=length(dirichlet);
ind=1;
for i=1:nnode
    const(ind)=dirichlet(i)*ndimf-1;
    ind=ind+1;
    const(ind)=dirichlet(i)*ndimf;
    ind=ind+1;
end

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
s.material = material;
lhs = LHSintegrator.create(s);
k=lhs.compute();
end