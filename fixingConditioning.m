% Mesh
mesh = UnitTriangleMesh(3,3);
mesh.plot()

% Level Set
nNodes = size(mesh.coord,1);
lvlSet = -1*ones(nNodes,1);
lvlSet(5,1) = 1;

% Unfitted Mesh
sUm.boundaryMesh   = mesh.createBoundaryMesh();
sUm.backgroundMesh = mesh;
uMesh = UnfittedMesh(sUm);
uMesh.compute(lvlSet);

figure()
uMesh.plot()

% Full inner mesh
full_inner = uMesh.createInnerMesh();
canonical = full_inner.computeCanonicalMesh();
figure()
canonical.plot