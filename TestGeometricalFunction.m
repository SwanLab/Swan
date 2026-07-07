clc,clear,close all

%% Load base mesh
% file = 'Tetradecahedron0025';%'Tetradecahedron0025';
% 
% TMC = TetradecahedronMeshComputer(file);
% mesh = TMC.getMesh();
% MS = TMC.getMasterSlave();

mesh = UnitTetraMesh(100,100,100)

%% Create level set
gPar.type = 'Tetradecahedron'; %'Octahedron';
gPar.xCoorCenter = 0.5;
gPar.yCoorCenter = 0.5;
gPar.zCoorCenter = 0.5;
gPar.radius = 0.25;

g         = GeometricalFunction(gPar);
phiFun    = g.computeLevelSetFunction(mesh);
lsCircle  = phiFun.fValues;
ls = -lsCircle;

%% Create unfitted mesh / density
sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh;
uMesh              = UnfittedMesh(sUm);
uMesh.compute(ls);
uMesh.plot

ls = CharacteristicFunction.create(uMesh);
s.trial = LagrangianFunction.create(mesh,1,'P1');
s.mesh = mesh;
f = FilterLump(s);
density = f.compute(ls,2);

%% set material
s.interpolation  = 'SIMPALL';
s.dim            = mesh.ndim;
s.matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(1e-5,0.3,mesh.ndim);
s.matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(1e-5,0.3);
s.matB.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(1,0.3,mesh.ndim);
s.matB.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(1,0.3);
mI = MaterialInterpolator.create(s);

x{1} = density;
s.mesh                 = mesh;
s.type                 = 'DensityBased';
s.density              = x;
s.materialInterpolator = mI;
s.dim = '3D';
material = Material.create(s);

%% set boundary conditions
isV1X = @(coor) (abs(coor(:,1) - 0.25) < 1e-5);
isV1Y = @(coor) (abs(coor(:,2) - 0.5) < 1e-5);
isV1Z = @(coor) (abs(coor(:,3) - 0) < 1e-5);
isVertex1 = @(coor) isV1X(coor) & isV1Y(coor) & isV1Z(coor);

sDir{1}.domain    = @(coor) isVertex1(coor);
sDir{1}.direction = [1,2,3];
sDir{1}.value     = 0;

isV2X = @(coor) (abs(coor(:,1) - 0.75) < 1e-5);
isV2Y = @(coor) (abs(coor(:,2) - 0.5) < 1e-5);
isV2Z = @(coor) (abs(coor(:,3) - 0) < 1e-5);
isVertex2 = @(coor) isV2X(coor) & isV2Y(coor) & isV2Z(coor);

sDir{2}.domain    = @(coor) isVertex2(coor);
sDir{2}.direction = [2];
sDir{2}.value     = 0;

isV3X = @(coor) (abs(coor(:,1) - 0.5) < 1e-5);
isV3Y = @(coor) (abs(coor(:,2) - 0.25) < 1e-5);
isV3Z = @(coor) (abs(coor(:,3) - 0) < 1e-5);
isVertex3 = @(coor) isV3X(coor) & isV3Y(coor) & isV3Z(coor);

sDir{3}.domain    = @(coor) isVertex3(coor);
sDir{3}.direction = [3];
sDir{3}.value     = 0;

dirichletFun = [];
for i = 1:numel(sDir)
    dir = DirichletCondition(mesh, sDir{i});
    dirichletFun = [dirichletFun, dir];
end
s.dirichletFun = dirichletFun;
s.pointloadFun = [];
s.periodicFun  = 1; %Set to not be empty
s.mesh = mesh;
bc = BoundaryConditions(s);
bc.updatePeriodicConditions(MS);

%% Set micro problem
s.mesh     = mesh;
s.material = material;
s.scale    = 'MICRO';
s.dim      = '3D';
s.boundaryConditions = bc;

%% Set solver
BCAp = BCApplier(s);

rigModes = RigidBodyFunction.create(mesh,[0.5,0.5,0.5]);
RFun = rigModes.projectBasisFunctions('P1');
for i = 1:length(RFun)
    Rfull(:,i) = reshape(RFun{i}.fValues',[],1);
end

for i = 1:size(Rfull,2)
    R(:,i) = BCAp.fullToReducedVectorDirichlet(Rfull(:,i));
end
s.type = 'ELASTIC';
s.nullSpace = R;
s.nLevels = 5;
s.tol = 1e-8;
s.maxIter = 1;
p     = pyAMG.create(s);

sS.preconditioner = p;
sS.tol = 1e-5;
solver = PCG(sS);

%% Continue problem definition

s.solverCase = solver;%DirectSolver();
s.solverType = 'REDUCED';
s.solverMode = 'FLUC';
fem = ElasticProblemMicro(s);
material.setDesignVariable({density})
fem.updateMaterial(material.obtainTensor())
fem.solve();

totVol = mesh.computeVolume();
matHomog = fem.Chomog/totVol;
phi = 1 - Integrator.compute(density,mesh,1)/(1/2);

C11 = matHomog(1,1,1,1);
C12 = matHomog(1,1,2,2);
C44 = matHomog(2,3,2,3);
ZenerRatioTetra = 2*C44./(C11-C12);
