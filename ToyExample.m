

clear;
close all;
clc;


% MESH
m = TriangleMesh(1,1,50,50);


% FIBER UNFITTED FUNCTION
sG.type        = 'HorizontalFiber';
sG.width       = 0.25;
sG.yCoorCenter = 0.5;
g              = GeometricalFunction(sG);
ls             = g.computeLevelSetFunction(m);

sU.backgroundMesh = m;
sU.boundaryMesh   = m.createBoundaryMesh();
uM                = UnfittedMesh(sU);
uM.compute(ls.fValues);

chi = CharacteristicFunction.create(uM);


% NONLINEAR SEGMENT
s.mesh  = m;
s.theta = 90;
s.alpha = 8;
s.beta  = 0;


% EPSILON STUDY
h = m.computeMeanCellSize();
epsVec = 1.5*h:0.5*h:5*h;
LineSearchVec = [0.01, 0.1, 1, 10];

for i = 1:length(LineSearchVec)
    s.lineSearch = LineSearchVec(i);
    filter  = NonLinearFilterSegment(s);
    for j = 1:length(epsVec)
        filter.updateEpsilon(epsVec(j));
        rhoEps{i,j} = filter.compute(chi,2);
        iter(i,j) = filter.iter;
    end
end



