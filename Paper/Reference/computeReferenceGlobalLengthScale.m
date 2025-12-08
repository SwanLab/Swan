function ls = computeReferenceGlobalLengthScale(matFile)
    load(matFile,'d');
    mesh = d.fun.mesh;
    h    = mesh.computeMeanCellSize();
    
    s.mesh = mesh;
    s.filterType = 'PDE';
    s.trial = LagrangianFunction.create(mesh,1,'P1');
    filter  = Filter.create(s);
    
    levelSet          = -ones(mesh.nnodes,1);
    sG.backgroundMesh = mesh;
    sG.boundaryMesh   = mesh.createBoundaryMesh();
    uMesh = UnfittedMesh(sG);
    uMesh.compute(levelSet);
    
    ss.mesh       = mesh;
    ss.uMesh      = uMesh;
    ss.filter     = filter;
    ss.minEpsilon = 3*h;
    ss.epsilon    = 3*h;
    ss.target     = 1;
    ss.value0     = 1;
    ss.test       = s.trial;
    lsF           = LengthScaleConstraint(ss);
    
    ls = lsF.computeFunctionAndGradient(d);
    ls = (-ls+1)*100;
end