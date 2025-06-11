fHandle = @(x) [x(1,:,:); x(1,:,:)+x(2,:,:); 2*x(2,:,:)];
mesh = QuadMesh(1,1,2,2);
aF = AnalyticalFunction2.create(fHandle,mesh);
aF.evaluate([0 1; 0 1])