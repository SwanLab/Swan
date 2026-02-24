clear;
close all;

x1       = linspace(0,1,30);
x2       = linspace(0,1,30);
[xv,yv]  = meshgrid(x1,x2);
[F,V]    = mesh2tri(xv,yv,zeros(size(xv)),'x');
m.coord  = V(:,1:2);
m.connec = F;
mesh     = Mesh.create(m);

gPar.type         = 'CrossedSquare';
gPar.length       = 1;
gPar.xCoorCenter  = 0.5;
gPar.yCoorCenter  = 0.5;
gPar.tFrame       = 0.01;
gPar.tCross       = 0.01;
% 
% gPar.type       = 'Square';
% gPar.length     =   0.5 ;
% gPar.xCoorCenter  =  0.5;
% gPar.yCoorCenter  =  0.5;

% gPar.type         = 'Circle';
% gPar.radius       = 0.25;
% gPar.xCoorCenter  = 0.5;
% gPar.yCoorCenter  = 0.5;

g                 = GeometricalFunction(gPar);

phiFun            = g.computeLevelSetFunction(mesh);
lsCircle          = phiFun.fValues;
lsCircleInclusion = -lsCircle;

sUm.backgroundMesh = mesh;
sUm.boundaryMesh   = mesh.createBoundaryMesh;
uMesh              = UnfittedMesh(sUm);
uMesh.compute(lsCircleInclusion);



            % 
            % [ls,phiFun] = obj.computeLevelSet(obj.mesh);
            % sUm.backgroundMesh = obj.mesh;
            % sUm.boundaryMesh   = obj.mesh.createBoundaryMesh;
            % uMesh              = UnfittedMesh(sUm);
            % uMesh.compute(ls);
            % %holeMesh = uMesh.createInnerMesh();
            % %obj.mesh = holeMesh;            
            % close all;
            % uMesh.plot;
            % funLS     = CharacteristicFunction.create(uMesh);          
            % s.filterType = 'LUMP';
            % s.mesh  = obj.mesh;
            % s.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            % f = Filter.create(s); 
            % obj.filter = f;
            % %obj.density = f.compute(funLS,2);
            % s.fun = f.compute(funLS,2);
            % s.type = 'Density';
            % s.plotting = true;
            % dens               = DesignVariable.create(s);
            % obj.designVariable = dens;
            % 




uMesh.plot;