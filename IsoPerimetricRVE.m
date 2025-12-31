% Isoperimetric RVE

%% Intro

clear;
clc;

epsOverh = 3;
mesh = TriangleMesh(1,1,100,100);
tau = 0.1;

V = createVolume(mesh);
P = createPerimeter(mesh,epsOverh);

%% Circle 1

s.type = 'Circle';
s.xCoorCenter = 0.5;
s.yCoorCenter = 0.5;
sD.mesh = mesh;
sD.type = 'LevelSet';
sD.plotting = false;
R = linspace(0,sqrt(2)/2,100);
for i = 1:length(R)
    s.radius = R(i);
    g        = GeometricalFunction(s);
    lsF      = g.computeLevelSetFunction(mesh);
    sD.fun   = lsF+1e-4;
    ls       = DesignVariable.create(sD);
    Vx       = V.computeFunctionAndGradient(ls);
    Px       = P.computeFunctionAndGradient(ls);
    rhoHC1(i) = Vx;
    IPC1(i)   = Vx/(Px+tau);
end

%% Circle 2

s.type = 'Circle';
s.xCoorCenter = 1.0;
s.yCoorCenter = 0.5;
sD.mesh = mesh;
sD.type = 'LevelSet';
sD.plotting = false;
R = linspace(0,sqrt(5)/2,100);
for i = 1:length(R)
    s.radius = R(i);
    g        = GeometricalFunction(s);
    lsF      = g.computeLevelSetFunction(mesh);
    sD.fun   = lsF+1e-4;
    ls       = DesignVariable.create(sD);
    Vx       = V.computeFunctionAndGradient(ls);
    Px       = P.computeFunctionAndGradient(ls);
    rhoHC2(i) = Vx;
    IPC2(i)   = Vx/(Px+tau);
end

%% Circle 3

s.type = 'Circle';
s.xCoorCenter = 1.0;
s.yCoorCenter = 0.0;
sD.mesh = mesh;
sD.type = 'LevelSet';
sD.plotting = false;
R = linspace(0,sqrt(2),100);
for i = 1:length(R)
    s.radius = R(i);
    g        = GeometricalFunction(s);
    lsF      = g.computeLevelSetFunction(mesh);
    sD.fun   = lsF+1e-4;
    ls       = DesignVariable.create(sD);
    Vx       = V.computeFunctionAndGradient(ls);
    Px       = P.computeFunctionAndGradient(ls);
    rhoHC3(i) = Vx;
    IPC3(i)   = Vx/(Px+tau);
end

%% Fiber 1

s.type = 'HorizontalFiber';
s.yCoorCenter = 0.5;
sD.mesh = mesh;
sD.type = 'LevelSet';
sD.plotting = false;
h = linspace(0,1,100);
for i = 1:length(R)
    s.width = h(i);
    g        = GeometricalFunction(s);
    lsF      = g.computeLevelSetFunction(mesh);
    sD.fun   = lsF+1e-4;
    ls       = DesignVariable.create(sD);
    Vx       = V.computeFunctionAndGradient(ls);
    Px       = P.computeFunctionAndGradient(ls);
    rhoHF1(i) = Vx;
    IPF1(i)   = Vx/(Px+tau);
end

%% Fiber 2

s.type = 'DiagonalFiber';
s.xVertex = 0.5;
s.yVertex = 0.0;
s.alpha = 0;
sD.mesh = mesh;
sD.type = 'LevelSet';
sD.plotting = false;
h = linspace(0,1,100);
for i = 1:length(R)
    s.width = h(i);
    g        = GeometricalFunction(s);
    lsF      = g.computeLevelSetFunction(mesh);
    sD.fun   = lsF+1e-4;
    ls       = DesignVariable.create(sD);
    Vx       = V.computeFunctionAndGradient(ls);
    Px       = P.computeFunctionAndGradient(ls);
    rhoHF2(i) = Vx;
    IPF2(i)   = Vx/(Px+tau);
end

%% Fiber 3

s.type = 'DiagonalFiber';
s.xVertex = 1.0;
s.yVertex = 0.0;
s.alpha = pi/4;
sD.mesh = mesh;
sD.type = 'LevelSet';
sD.plotting = false;
h = linspace(0,sqrt(2),100);
for i = 1:length(R)
    s.width = h(i);
    g        = GeometricalFunction(s);
    lsF      = g.computeLevelSetFunction(mesh);
    sD.fun   = lsF+1e-4;
    ls       = DesignVariable.create(sD);
    Vx       = V.computeFunctionAndGradient(ls);
    Px       = P.computeFunctionAndGradient(ls);
    rhoHF3(i) = Vx;
    IPF3(i)   = Vx/(Px+tau);
end

%% N Fibers

s.type = 'HorizontalNFibers';
s.nFibers = 4;
s.minyCoor = 0;
s.maxyCoor = 1;
sD.mesh = mesh;
sD.type = 'LevelSet';
sD.plotting = false;
h = linspace(0,1/4,100);
for i = 1:length(R)
    s.width = h(i);
    g        = GeometricalFunction(s);
    lsF      = g.computeLevelSetFunction(mesh);
    sD.fun   = lsF+1e-4;
    ls       = DesignVariable.create(sD);
    Vx       = V.computeFunctionAndGradient(ls);
    Px       = P.computeFunctionAndGradient(ls);
    rhoHFs(i) = Vx;
    IPFs(i)   = Vx/(Px+tau);
end

%% Functions

function V = createVolume(mesh)
    levelSet         = -ones(mesh.nnodes,1);
    s.backgroundMesh = mesh;
    s.boundaryMesh   = mesh.createBoundaryMesh();
    uMesh = UnfittedMesh(s);
    uMesh.compute(levelSet);

    s.mesh   = mesh;
    s.test = LagrangianFunction.create(mesh,1,'P1');
    s.uMesh = uMesh;
    V = VolumeFunctional(s);
end

function P = createPerimeter(mesh,epsOverh)
    levelSet         = -ones(mesh.nnodes,1);
    s.backgroundMesh = mesh;
    s.boundaryMesh   = mesh.createBoundaryMesh();
    uMesh = UnfittedMesh(s);
    uMesh.compute(levelSet);

    s.mesh  = mesh;
    s.trial = LagrangianFunction.create(mesh,1,'P1');
    s.filterType  = 'PDE';
    filter = Filter.create(s);

    h         = mesh.computeMeanCellSize();
    s.mesh    = mesh;
    s.uMesh   = uMesh;
    s.filter  = filter;
    s.epsilon = epsOverh*h;
    s.value0  = 1;
    P         = PerimeterFunctional(s);
end