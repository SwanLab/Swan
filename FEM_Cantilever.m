%function FEM_Cantilever

clc;clear;close all

addpath(genpath(fileparts(mfilename('fullpath'))))

%file = 'test2d_triangle';
%a.fileName = file;
% s = FemDataContainer(a);


% Generate coordinates
x1 = linspace(0,7,15);
x2 = linspace(1,2,15);
% Create the grid
[xv,yv] = meshgrid(x1,x2);
% Triangulate the mesh to obtain coordinates and connectivities
[F,V] = mesh2tri(xv,yv,zeros(size(xv)),'x');

s.coord = V(:,1:2);
s.connec = F;
mesh = Mesh(s);

s.mesh = mesh;
s.type = 'ELASTIC';
s.scale = 'MACRO';
s.material = createMaterial(mesh,1);
s.dim = '2D';
s.bc = createBoundaryConditions(mesh);
fem = FEM.create(s);
fem.solve();

figure
fem.uFun.plot()
figure
fem.stressFun.plot()
figure
fem.strainFun.plot()

fem.uFun.fValues(:,end+1) = 0;
fem.uFun.ndimf = 3;

%fem.uFun.print('results_fem_disp', 'Paraview') % print using Paraview
%fem.print('results_fem2', 'Paraview') % print using Paraview

%plotError()

%end

function bc = createBoundaryConditions(mesh)
%     dirichletNodes = abs(mesh.coord(:,1)-0) < 1e-12;
%     rightSide  = max(mesh.coord(:,1));
%     isInRight = abs(mesh.coord(:,1)-rightSide)< 1e-12;
%     isInMiddleEdge = abs(mesh.coord(:,2)-1.5) < 0.1;
%     forceNodes = isInRight & isInMiddleEdge;
%     nodes = 1:mesh.nnodes;
%     bcDir = [nodes(dirichletNodes)';nodes(dirichletNodes)'];
%     bcDir(1:2:end,end+1) = 1;
%     bcDir(2:2:end,end) = 2;
%     bcDir(:,end+1) = 0;
%     bc.dirichlet = bcDir;
%     bc.dirichlet = nodes(dirichletNodes);
%     bc.pointload(:,1) = nodes(forceNodes);
%     bc.pointload(:,2) = 2;
%     bc.pointload(:,3) = -1;
%     
    dirichletNodes = abs(mesh.coord(:,1)-0) < 1e-12;
    % rightSide  = max(mesh.coord(:,1));
    % isInRight = abs(mesh.coord(:,1)-rightSide)< 1e-12;
    % isInMiddleEdge = abs(mesh.coord(:,2)-1.5) < 0.1;
    topSide = max(mesh.coord(1,:));
    isOnTop = abs(mesh.coord(:,1)-topSide)< 1e-12;
    forceNodes = isOnTop;
    nodes = 1:mesh.nnodes;
    bcDir = [nodes(dirichletNodes)';nodes(dirichletNodes)'];
    nodesdir=size(nodes(dirichletNodes),2);
    bcDir(1:nodesdir,end+1) = 1;
    bcDir(nodesdir+1:end,end) = 2;
    bcDir(:,end+1)=0;
    bc.dirichlet = bcDir;
    bc.pointload(:,1) = nodes(forceNodes);
    bc.pointload(:,2) = 2;

    n = 15;
    x = linspace(0,13,n);
    y = zeros(1,length(x));
    suma = 0;
    for i = 1:length(x)
        y(i) = sqrt((1 - x(i)^2/(13^2)) * (0.52^2));
        suma = suma + y(i);
    end
    v =  150;
    S = 55/2;
    rho = 1.225 *((288.15-6.5e-3*6500)/288.15)^(-1+9.81/(287*6.5e-3));    
    y = 1/8 * rho * v^2 * S*y; %Lift
    for i = 1:length(nodes(forceNodes))-1
        bc.pointload(i,3) = -y(i);
    end
end

function material = createMaterial(mesh,ngaus)
    s.mesh = mesh;
    s.type = 'ELASTIC';
    s.scale = 'MACRO';
    ngaus = 1;
    I = ones(mesh.nelem,ngaus);
    s.ptype = 'ELASTIC';
    s.pdim  = '2D';
    s.nelem = mesh.nelem;
    s.mesh  = mesh;
    E = 70e9; 
    nu = 0.35;
    kappa = E/(2*(1-nu));
    s.kappa = kappa*I;
    %mu = (kappa - kappa*nu)/(1 + nu);
    mu = E./(2*(1+nu));
    s.mu    = mu*I;
    mat = Material.create(s);
    mat.compute(s);
    material = mat;
end