
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
epsilon = 4*h;
sF.trial = LagrangianFunction.create(mesh,1,'P1');
sF.mesh = mesh;
sF.filterType = 'PDE';
sF.metric = 'Isotropy';
filter = Filter.create(sF);
filter.updateEpsilon(epsilon);
rhoEps1 = filter.compute(chi,3);  % Cubic Gauss quadrature to integrate
rhoEps1.plot();


% Filter PDE Anisotropy



% Filter Lump



% (Optional) Is not in the tutorial !!   NonLinearFilterSegment

