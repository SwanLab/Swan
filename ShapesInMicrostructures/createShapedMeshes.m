function createShapedMeshes
% Define the initial data
D0.c = [0.6204,0.6204,0.6204];
D0.theta = [0,60,120];
D0.divUnit = 100;
D0.filename = 'HexagonTest.m';
MC = MeshCreator(D0);
MC.computeMeshNodes();
MC.drawMesh();
MC.plotIndicesOfNodes();
end