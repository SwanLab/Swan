load('negPoissAmitges.mat')
background = Mesh.createFromGiD('test2d_micro');
boundary = background.createBoundaryMesh();

a.backgroundMesh = background;
a.boundaryMesh   = boundary;
unfitted = UnfittedMesh(a);
unfitted.compute(levelsetintermig.fun.fValues)
IM = unfitted.createInnerMesh();