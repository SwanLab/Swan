function exportT(T,R,mesh,name)

z.mesh      = mesh;
z.order     = 'P1';

for i=1:8
  z.fValues   = reshape(T(:,i),[mesh.ndim,mesh.nnodes])';
  uFeFun = LagrangianFunction(z);%
  fileName = [name num2str(i)];
  centroids=computeCentroid(mesh);
  CoarsePlotSolution(uFeFun, mesh,[],fileName, R, centroids);
  %uFeFun.print(fileName,'Paraview');
end

end


function [centroids] =computeCentroid(mesh)
    x0=mean(mesh.coord(:,1));
    y0=mean(mesh.coord(:,2));
    centroids = [x0,y0];
end