function exportT_holeInclusion(T,mesh,name)
    z.mesh      = mesh;
    z.order     = 'P1';
    
    for i=1:size(T,2)
      z.fValues   = reshape(T(:,i),[mesh.ndim,mesh.nnodes])';
      uFeFun = LagrangianFunction(z);%
      fileName = name +num2str(i);
      uFeFun.print(fileName,'Paraview');
    end
end
