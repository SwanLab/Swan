% This script has the purpose to test the geometrical functions modified

%% Main
mesh=createMesh();
type= 'Circle';

[ls,phiFun] = computeLevelSet(mesh,type);
sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh;
uMesh              = UnfittedMesh(sUm);
uMesh.compute(lsCircleInclusion);
uMesh.plot;



%% Functions

function mesh=createMesh()
    mesh = UnitTriangleMesh(20,20);
    % mesh = UnitTetraMesh(20,20,20);
end

function [ls,phiFun] = computeLevelSet(mesh,type)            
     g.type          = type;
     g.minxCoor      = 0;
     g.maxxCoor      = 1;
     g.minyCoor      = 0;
     g.maxyCoor      = 1;
     g.minzCoor      = 0;
     g.maxzCoor      = 1; 
     g.width         = 0.5;
     g.nFibersY      = 4;    
     g.nFibersZ      = 4;
     g.radius        = 0.075;
     g               = GeometricalFunction(g);
     phiFun          = g.computeLevelSetFunction(mesh);
     lsValues        = phiFun.fValues;
     ls              = lsValues;
end