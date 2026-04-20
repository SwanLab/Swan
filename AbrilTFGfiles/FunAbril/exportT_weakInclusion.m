function exportT_weakInclusion(T,R,mesh,uMesh,name)

z.mesh      = mesh;
z.order     = 'P1';

for i=1:size(T,2)
  z.fValues   = reshape(T(:,i),[mesh.ndim,mesh.nnodes])';
  uFeFun = LagrangianFunction(z);%
  
  uMeshFun = uMesh.obtainFunctionAtUnfittedMesh(uFeFun);
  fvalues = [uMeshFun.innerMeshFunction.fValues;
      uMeshFun.innerCutMeshFunction.fValues];

  s.coord = [uMeshFun.innerMeshFunction.mesh.coord;
      uMeshFun.innerCutMeshFunction.mesh.coord];

  s.connec = [uMeshFun.innerMeshFunction.mesh.connec;
      uMeshFun.innerCutMeshFunction.mesh.connec  + max(uMeshFun.innerMeshFunction.mesh.connec(:))];

  mh = Mesh.create(s);
  ss.mesh = mh;
  ss.fValues = fvalues;
  ss.order = 'P1';
  ss.ndimf = size(fvalues,2)
  u = LagrangianFunction(ss);

  fileName = strrep("r" + num2str(R), '.', '_')+ name +num2str(i);
  u.print(fileName,'Paraview')

end

end


function [centroids] =computeCentroid(mesh)
centroids = zeros(size(mesh.coord,2),1);
for i = 1:size(mesh.coord,2)
    centroids(i,1) =mean(mesh.coord(:,i));

end
end