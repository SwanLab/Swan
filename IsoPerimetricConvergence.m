% Isoperimetric length scale convergence

clear;
clc;

matFile = 'Paper/Reference/DensityVerticalCantilever.mat';
load(matFile,'d');
mesh = d.fun.mesh;
h    = mesh.computeMeanCellSize();

s.mesh = mesh;
s.filterType = 'PDE';
s.trial = LagrangianFunction.create(mesh,1,'P1');
filter  = Filter.create(s);

ss.mesh       = mesh;
ss.filter     = filter;
ss.minEpsilon = 3*h;
ss.epsilon    = 3*h;
ss.target     = 0.05;
ss.value0     = 1;
ss.test       = s.trial;

K = 2.^((1:5)-1);
for k = 1:length(K)
    p0Mesh = QuadMesh(1,2,K(k),K(k));
    lsFun{k} = LagrangianFunction.create(p0Mesh,1,'P0');
    N = K(k);
    lsVal = [];
    for i = 1:N
        for j = 1:N
            x0       = i/N - 1/(2*N);
            y0       = 2*j/N - 1/N;
            ss.uMesh = createBaseDomain(mesh,x0,y0,N);
            ss.delta = 0.1/N;
            lsF      = LengthScaleConstraint(ss);
            ls       = lsF.computeFunctionAndGradient(d);
            ls       = (-ls+1)*100*0.05;
            lsVal    = [lsVal;ls];
        end
    end
    lsFun{k}.setFValues(double(lsVal<0.05));
end






function uMesh = createBaseDomain(mesh,x0,y0,N)
    s.type             = 'Rectangle';
    s.xCoorCenter      = x0;
    s.yCoorCenter      = y0;
    s.xSide            = 1/N;
    s.ySide            = 2/N;
    g                  = GeometricalFunction(s);
    lsFun              = g.computeLevelSetFunction(mesh);
    sUm.backgroundMesh = mesh;
    sUm.boundaryMesh   = mesh.createBoundaryMesh();
    uMesh              = UnfittedMesh(sUm);
    uMesh.compute(lsFun.fValues);
end