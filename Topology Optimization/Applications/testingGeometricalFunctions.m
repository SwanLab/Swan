clear
close all

% Circle

file = 'Circle_Triangle_Linear_Unstructured';
a.fileName = file;
ff = FemDataContainer(a);
mesh = ff.mesh;
bMesh = mesh.createBoundaryMesh();

coor = mesh.coord;

s.type='Circle';
s.radius=0.3*(max(coor(:,1))-min(coor(:,1)));
s.xCoorCenter=0.5*(max(coor(:,1))-min(coor(:,1)));
s.yCoorCenter=0.5*(max(coor(:,2))-min(coor(:,2)));
g=GeometricalFunction(s);

ls=g.computeLevelSetFunction(mesh);

ss.backgroundMesh = mesh;
ss.boundaryMesh   = bMesh;
uMesh             = UnfittedMesh(ss);
uMesh.compute(ls.fValues);
uMesh.plot();




clear

% Square

file = 'Circle_Triangle_Linear_Unstructured';
a.fileName = file;
ff = FemDataContainer(a);
mesh = ff.mesh;
bMesh = mesh.createBoundaryMesh();

coor = mesh.coord;

s.type='Square';
s.length=0.6*(max(coor(:,1))-min(coor(:,1)));
s.xCoorCenter=0.5*(max(coor(:,1))-min(coor(:,1)));
s.yCoorCenter=0.5*(max(coor(:,2))-min(coor(:,2)));
g=GeometricalFunction(s);

ls=g.computeLevelSetFunction(mesh);

ss.backgroundMesh = mesh;
ss.boundaryMesh   = bMesh;
uMesh             = UnfittedMesh(ss);
uMesh.compute(ls.fValues);
figure
uMesh.plot();



clear

% Fiber

file = 'Circle_Triangle_Linear_Unstructured';
a.fileName = file;
ff = FemDataContainer(a);
mesh = ff.mesh;
bMesh = mesh.createBoundaryMesh();

coor = mesh.coord;

s.type='VerticalFiber';
s.width=0.3*(max(coor(:,1))-min(coor(:,1)));
s.xCoorCenter=0.5*(max(coor(:,1))-min(coor(:,1)));
g=GeometricalFunction(s);

ls=g.computeLevelSetFunction(mesh);

ss.backgroundMesh = mesh;
ss.boundaryMesh   = bMesh;
uMesh             = UnfittedMesh(ss);
uMesh.compute(ls.fValues);
figure
uMesh.plot();




clear

% Vertical N-Fibers

file = 'Circle_Triangle_Linear_Unstructured';
a.fileName = file;
ff = FemDataContainer(a);
mesh = ff.mesh;
bMesh = mesh.createBoundaryMesh();

coor = mesh.coord;

s.type='VerticalNFibers';
s.nFibers=5;
s.minxCoor=min(coor(:,1));
s.maxxCoor=max(coor(:,1));
g=GeometricalFunction(s);

ls=g.computeLevelSetFunction(mesh);

ss.backgroundMesh = mesh;
ss.boundaryMesh   = bMesh;
uMesh             = UnfittedMesh(ss);
uMesh.compute(ls.fValues);
figure
uMesh.plot();