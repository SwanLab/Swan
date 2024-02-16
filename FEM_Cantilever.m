function FEM_Cantilever

clc;clear;close all

addpath(genpath(fileparts(mfilename('fullpath'))))

%file = 'test2d_triangle';
%a.fileName = file;
% s = FemDataContainer(a);


% Generate coordinates
x1 = linspace(0,2,20);
x2 = linspace(1,2,20);
% Create the grid
[xv,yv] = meshgrid(x1,x2);
% Triangulate the mesh to obtain coordinates and connectivities
[F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');

s.coord = V(:,1:2);
s.connec = F;
mesh = Mesh.create(s);

s.mesh = mesh;
s.type = 'ELASTIC';
s.scale = 'MACRO';
s.material = createMaterial(mesh,1);
s.dim = '2D';
s.bc = createBoundaryConditions(mesh);
fem = FEM.create(s);
fem.solve();

figure(1)
fem.uFun.plot()
figure(2)
fem.stressFun.plot()
figure(3)
fem.strainFun.plot()

fem.print('results_fem', 'Paraview') % print using Paraview

end

function bc = createBoundaryConditions(mesh)
dirichletNodes = abs(mesh.coord(:,1)-0) < 1e-12;
rightSide  = max(mesh.coord(:,1));
isInRight = abs(mesh.coord(:,1)-rightSide)< 1e-12;
isInMiddleEdge = abs(mesh.coord(:,2)-1.5) < 0.1;
forceNodes = isInRight & isInMiddleEdge;
nodes = 1:mesh.nnodes;
bc.dirichlet = nodes(dirichletNodes);
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