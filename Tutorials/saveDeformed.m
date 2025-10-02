function deformed = saveDeformed(mesh,fvalues)
    
    s.mesh  = mesh;
    s.ndimf = 2;
    s.order = 'P1';
    s.fValues = reshape(fvalues,2,[])' ;
    deformed   = LagrangianFunction(s);
end