clc,clear,close all

%% Load base mesh
file = 'Tetradecahedron';

TMC = TetradecahedronMeshComputer(file);
mesh = TMC.getMesh();
MS = TMC.getMasterSlave();

%% Create level set
gPar.type = 'Tetradecahedron'; %'Octahedron';
gPar.xCoorCenter = 0.5;
gPar.yCoorCenter = 0.5;
gPar.zCoorCenter = 0.5;
gPar.radius = 0.75;

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
s.matA.bulk  = IsotropicElasticMaterial.computeKappaFromYoungAndPoisson(1e-12,0.3,mesh.ndim);
s.matA.shear = IsotropicElasticMaterial.computeMuFromYoungAndPoisson(1e-12,0.3);
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
isCenterX = @(coor) (abs(coor(:,1) - 0.75) < 1e-5);
isCenterY = @(coor) (abs(coor(:,2) - 0.5) < 1e-5);
isCenterZ = @(coor) (abs(coor(:,3) - 0) < 1e-5);
sDir{1}.domain    = @(coor) isCenterX(coor) & isCenterY(coor) & isCenterZ(coor);
sDir{1}.direction = [1,2,3];
sDir{1}.value     = 0;

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
%bc.updatePeriodicConditions(MS);

%% Set micro problem
s.mesh     = mesh;
s.material = material;
s.scale    = 'MICRO';
s.dim      = '3D';
s.boundaryConditions = bc;
s.solverCase = DirectSolver();
s.solverType = 'REDUCED';
s.solverMode = 'FLUC';
fem = ElasticProblemMicro(s);
material.setDesignVariable({density})
fem.updateMaterial(material.obtainTensor())
fem.solve();

totVol = mesh.computeVolume();
matHomog = fem.Chomog/totVol;
phi = 1 - Integrator.compute(obj.density,obj.baseMesh,1)/(1/2);

C11 = matHomog(1,1,1,1);
C12 = matHomog(1,1,2,2);
C44 = matHomog(2,3,2,3);
ZenerRatioTetra = 2*C44./(C11-C12);
