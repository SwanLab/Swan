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
x1 = linspace(0,1,5);
x2 = linspace(0,1,5);
% Create the grid
[xv,yv] = meshgrid(x1,x2);
% Triangulate the mesh to obtain coordinates and connectivities
[F,basis] = mesh2tri(xv,yv,zeros(size(xv)),'x');

s.coord = basis(:,1:2);
s.connec = F;
mesh = Mesh(s);

s.mesh = mesh;
s.type = 'ELASTIC';
s.scale = 'MACRO';
s.material = createMaterial(mesh,1);
s.dim = '2D';
s.bc = createBoundaryConditions(mesh);
nDimf=2;

dof=1:1:nDimf*mesh.nnodes;
const = constDof(s.bc.dirichlet,nDimf);
freedof=setdiff(dof,const);

dispFun=P1Function.create(s.mesh, nDimf);
K= computeStiffnessMatrix(s.mesh,s.material,dispFun);
Kred=K(freedof,freedof);
[basis,D]=eig(full(Kred));

fvalues=zeros(nDimf*mesh.nnodes,1);
fvalues(s.bc.dirichlet)=0;
fvalues(freedof)=basis(:,1);
fvalues=reshape(fvalues,[mesh.ndim,mesh.nnodes])';

dispFun.fValues=fvalues;
dispFun.plot
fem = FEM.create(s);
fem.solve();

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








function bc = createBoundaryConditions(mesh)
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