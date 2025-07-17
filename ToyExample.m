

clear;
close all;
clc;


% MESH
m = TriangleMesh(1,1,150,150);


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

% fValues = 1-heaviside(ls.fValues);
% chi = LagrangianFunction.create(m,1,'P1');
% chi.setFValues(fValues);


% NONLINEAR SEGMENT
s.mesh  = m;
s.theta = 90;
s.alpha = 4;
s.beta  = 0;


% EPSILON STUDY
h = m.computeMeanCellSize();
epsVec = 2*h; % 1.5*h:0.5*h:5*h                [2*h,5*h,10*h,20*h]
%LineSearchVec = 500; % [0.01, 0.1, 1, 10, 100]           [0.01, 0.1, 1, 10]

%for i = 1:length(LineSearchVec)
    %s.lineSearch = LineSearchVec(i);
    for j = 1:length(epsVec)
        filter  = NonLinearFilterSegment(s);
        filter.updateEpsilon(epsVec(j));
        rhoEps{j} = filter.compute(chi,2);
        iter(j) = filter.iterVec(end);
        figure
        plot(filter.iterVec,filter.errorVec);
        %title(['LS',num2str(LineSearchVec(i)),'Eps',num2str(epsVec(j))]);
        title(['Eps',num2str(epsVec(j))]);
    end
%end



