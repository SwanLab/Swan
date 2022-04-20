function createShapedMeshes
% Define the initial data
D0.c = [2,6];
D0.theta = [0,102];
D0.divUnit = 1;
D0.filename = 'test.m';
MC = MeshCreator(D0);
MC.computeMeshNodes();
MC.drawMesh();
MC.plotIndicesOfNodes();
end