function plotRes(res,mesh,bc,malla,numero)
    xFull = bc.reducedToFullVector(res);
    s.fValues = reshape(xFull,2,[])';
    s.mesh = mesh;
    s.fValues(:,end+1) = 0;
    s.ndimf = 3;
    xF = P1Function(s);
    %xF.plot();
    xF.print(['Res',malla, num2str(numero)],'Paraview')
    fclose('all');
end