function createShapedMeshes
% Define the initial data
D0.c = [10/7,7/10];
D0.theta = [0,90];
D0.divUnit = 25;
D0.filename = 'RectangleTest.m';
MC = MeshCreator(D0);
MC.computeMeshNodes();
MC.drawMesh();
MC.plotIndicesOfNodes();
end