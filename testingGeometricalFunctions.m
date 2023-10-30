clear

file = 'test2d_triangle';
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