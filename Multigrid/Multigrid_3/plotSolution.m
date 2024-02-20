function plotSolution(x,mesh,bc, malla, numero)
    xFull = bc.reducedToFullVector(x);
    s.fValues = reshape(xFull,2,[])';
    s.mesh = mesh;
    s.fValues(:,end+1) = 0;
    s.ndimf = 3;
    xF = P1Function(s);
    xF.plot();
    %xF.print([malla,num2str(numero)],'Paraview')
    fclose('all');
end