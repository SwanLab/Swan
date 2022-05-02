function createShapedMeshes
% Define the initial data
D0.c = [1,1];
D0.theta = [45,135];
D0.divUnit = 25;
D0.filename = 'DiamondTest.m';
MC = MeshCreator(D0);
MC.computeMeshNodes();
MC.drawMesh();
MC.plotIndicesOfNodes();
end