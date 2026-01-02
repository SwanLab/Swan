function exportT_holeInclusion(T,R,mesh,name)
    z.mesh      = mesh;
    z.order     = 'P1';
    
    for i=1:8
      z.fValues   = reshape(T(:,i),[mesh.ndim,mesh.nnodes])';
      uFeFun = LagrangianFunction(z);%
      fileName = strrep("r" + num2str(R), '.', '_')+ name +num2str(i);
      uFeFun.print(fileName,'Paraview');
    end
    
end
