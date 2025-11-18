function P = computeReferenceGlobalPerimeterCircle(matFile)
    load(matFile,'d');
    mesh = d.fun.mesh;
    h    = mesh.computeMeanCellSize();
    
    s.mesh = mesh;
    s.filterType = 'PDE';
    s.trial = LagrangianFunction.create(mesh,1,'P1');
    filter = Filter.create(s);
    
    levelSet         = -ones(mesh.nnodes,1);
    sG.backgroundMesh = mesh;
    sG.boundaryMesh   = mesh.createBoundaryMesh();
    uMesh = UnfittedMesh(sG);
    uMesh.compute(levelSet);
    
    ss.mesh    = mesh;
    ss.uMesh   = uMesh;
    ss.filter  = filter;
    ss.epsilon = 2*h;
    ss.value0  = 1;
    pF         = PerimeterFunctional(ss);
    
    P          = pF.computeFunctionAndGradient(d);
end