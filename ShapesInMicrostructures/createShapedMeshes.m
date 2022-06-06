function createShapedMeshes
% Define the initial data
D0.c = [1,1];
D0.theta = [0,90];
D0.divUnit = 3;
D0.filename = 'SquareTest.m';
MC = MeshCreator(D0);
MC.computeMeshNodes();
MC.drawMesh();
MC.plotIndicesOfNodes();
end