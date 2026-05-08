
clear;
clc;

% Mesh

mesh = TriangleMesh(1,1,75,75);
h = mesh.computeMeanCellSize();


% Geometrical Fun

s.type = 'Circle';
s.radius = 0.3;
s.xCoorCenter = 0.5;
s.yCoorCenter = 0.5;
g = GeometricalFunction(s);
lsF = g.computeLevelSetFunction(mesh);


sUm.backgroundMesh = mesh;
sUm.boundaryMesh = mesh.createBoundaryMesh();
uMesh = UnfittedMesh(sUm);
uMesh.compute(lsF.fValues);
chi = CharacteristicFunction.create(uMesh);


% Filter PDE Isotropic
% epsilon = 4*h;
% sF.trial = LagrangianFunction.create(mesh,1,'P1');
% sF.mesh = mesh;
% sF.filterType = 'PDE';
% sF.metric = 'Isotropy';
% filter = Filter.create(sF);
% filter.updateEpsilon(epsilon);
% rhoEps1 = filter.compute(chi,3);  % Cubic Gauss quadrature to integrate
% rhoEps1.plot();


% Filter PDE Anisotropy

epsilon = 4*h;

sF = [];
sF.trial = LagrangianFunction.create(mesh,1,'P1');
sF.mesh = mesh;
sF.filterType = 'PDE';
sF.boundaryType = 'Neumann';
sF.metric = 'Anisotropy';
k = 10;
A = [k 0; 
     0  1/k];
R = [cosd(45), -sind(45); 
    sind(45), cosd(45)];

sF.A = ConstantFunction.create(R'*A*R,mesh);

filter = Filter.create(sF);
filter.updateEpsilon(epsilon);

rhoEps2 = filter.compute(chi,3);
rhoEps2.plot();

% Filter Lump

% sF = [];
% sF.trial = LagrangianFunction.create(mesh,1,'P1');
% sF.mesh = mesh;
% sF.filterType = 'LUMP';
% 
% filter = Filter.create(sF);
% 
% rhoLump = filter.compute(chi,3);
% rhoLump.plot();

% (Optional) Is not in the tutorial !!   NonLinearFilterSegment

% epsilon = 4*h;
% 
% sF = [];
% sF.filterType = 'Segment';
% sF.mesh       = mesh;
% 
% sF.alpha = 5;
% sF.beta  = 5;
% sF.theta = 45;   
% 
% filter = Filter.create(sF);
% filter.updateEpsilon(epsilon);
% 
% rhoSeg = filter.compute(chi,3);
% 
% figure;
% rhoSeg.plot();
% title('Nonlinear Segment Filter');
