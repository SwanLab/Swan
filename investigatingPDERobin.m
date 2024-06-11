clear;

% Mesh
m = HexaMesh(1,1,1,2,2,2);

% Density example
sD.fHandle = @(x) ones(size(x(1,:,:)));
sD.ndimf   = 1;
sD.mesh    = m;
aFun       = AnalyticalFunction(sD);
rho        = aFun.project('P1');

% Filter
sF.filterType   = 'PDE';
sF.mesh         = m;
sF.boundaryType = 'Robin';
sF.metric       = 'Isotropy';
sF.trial        = LagrangianFunction.create(m,1,'P1');
filter          = Filter.create(sF);

% Rho epsilon
rhoEps = filter.compute(rho,'LINEAR');

% Idea HEXA
% Given 1 Iso-hexa element with connec [1 2 3 4 5 6 7 8]:
% Quad 1: [1 2 3 4]
% Quad 2: [5 6 7 8]
% Quad 3: [1 2 6 5]
% Quad 4: [3 4 8 7]
% Quad 5: [2 3 7 6]
% Quad 6: [1 4 8 5]

bM1=m.createBoundaryMesh{1}.mesh;