clear;
clc;
close all;

mesh = createMesh();

E     = 1;
nu    = 0.49;
EFun  = AnalyticalFunction.create(@(x) E*ones(size(squeeze(x(1,:,:)))),1,mesh);
nuFun = AnalyticalFunction.create(@(x) nu*ones(size(squeeze(x(1,:,:)))),1,mesh);
EFun  = EFun.project('P1');
nuFun = nuFun.project('P1');

s         = [];
s.type    = 'ISOTROPIC';
s.ptype   = 'ELASTIC';
s.ndim    = mesh.ndim;
s.young   = EFun;
s.poisson = nuFun;
C         = Material.create(s);

s                    = [];
s.mesh               = mesh;
s.scale              = 'MACRO';
s.material           = C;
s.dim                = '3D';
s.boundaryConditions = createBoundaryConditions(mesh);
s.solverType         = 'REDUCED';
s.solverMode         = 'DISP';
fem                  = ElasticProblem(s);
fem.solve();

N = 3;
u = fem.uFun;

s           = [];
s.operation = @(xV) EFun.evaluate(xV)./(2*(1+nuFun.evaluate(xV)));
mu          = DomainFunction(s);

s           = [];
s.operation = @(xV) EFun.evaluate(xV)./(N*(1-(N-1)*nuFun.evaluate(xV)));
kappa       = DomainFunction(s);

% Compliance
strain      = SymGrad(u);
stress      = DDP(C,strain);
dCompliance = DDP(strain,stress);
c           = Integrator.compute(dCompliance,mesh,'QUADRATIC');

% Bulk compliance:
divu        = Divergence(u);
dbC         = DDP(kappa,DDP(divu,divu));
s           = [];
s.operation = @(xV) dbC.evaluate(xV)./c;
dbC         = DomainFunction(s);
bC          = Integrator.compute(dbC,mesh,'QUADRATIC');

% Shear compliance
e           = AntiVoigt(strain);
D           = Voigt(Deviatoric(e));
A           = VoigtDeviatorNormMaterial(mesh);
dsC         = DDP(mu,DDP(D,DDP(A,D)));
s           = [];
s.operation = @(xV) dsC.evaluate(xV)./c;
dsC         = DomainFunction(s);
sC          = Integrator.compute(dsC,mesh,'QUADRATIC');

% Bulk-Shear
NRG = dbC-dsC;
NRG = NRG.project('P1',mesh);

% Outputs:
disp(['Bulk: ',char(string(bC))]);
disp(['Shear: ',char(string(sC))]);
NRG.print('SIMPALL3DPaper/BenchmarkAnalysis/CantileverEnergyDistribution')










% THINGS TO EDIT:

function mesh = createMesh()
    mesh = HexaMesh(2,1,1,50,25,25);
end

function bc = createBoundaryConditions(mesh)
    xMax    = max(mesh.coord(:,1));
    yMax    = max(mesh.coord(:,2));
    zMax    = max(mesh.coord(:,3));
    isDir   = @(coor)  coor(:,1)==0;
    isForce = @(coor)  coor(:,1)==xMax & coor(:,2)>=0.3*yMax & ...
        coor(:,2)<=0.7*yMax & coor(:,3)>=0.3*zMax & coor(:,3)<=0.7*zMax;
    
    sDir{1}.domain    = @(coor) isDir(coor);
    sDir{1}.direction = [1,2,3];
    sDir{1}.value     = 0;
    
    sPL{1}.domain    = @(coor) isForce(coor);
    sPL{1}.direction = 3;
    sPL{1}.value     = -1;
    
    dirichletFun = [];
    for i = 1:numel(sDir)
        dir = DirichletCondition(mesh, sDir{i});
        dirichletFun = [dirichletFun, dir];
    end
    s.dirichletFun = dirichletFun;
    
    pointloadFun = [];
    for i = 1:numel(sPL)
        pl = PointLoad(mesh, sPL{i});
        pointloadFun = [pointloadFun, pl];
    end
    s.pointloadFun = pointloadFun;
    
    s.periodicFun  = [];
    s.mesh = mesh;
    bc = BoundaryConditions(s);
end