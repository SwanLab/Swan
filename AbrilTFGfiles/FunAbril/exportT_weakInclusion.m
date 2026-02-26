function exportT_weakInclusion(T,R,mesh,name)

z.mesh      = mesh;
z.order     = 'P1';

for i=1:size(T,2)
  z.fValues   = reshape(T(:,i),[mesh.ndim,mesh.nnodes])';
  uFeFun = LagrangianFunction(z);%
  fileName = strrep("r" + num2str(R), '.', '_')+ name +num2str(i);
  centroids=computeCentroid(mesh);
  CoarsePlotSolution(uFeFun, mesh,[],fileName, R, centroids);
  %uFeFun.print(fileName,'Paraview');
end

end


function [centroids] =computeCentroid(mesh)
centroids = zeros(size(mesh.coord,2),1);
for i = 1:size(mesh.coord,2)
    centroids(i,1) =mean(mesh.coord(:,i));

end
end