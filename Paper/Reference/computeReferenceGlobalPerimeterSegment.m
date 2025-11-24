function P = computeReferenceGlobalPerimeterSegment(matFile)
    load(matFile,'d');
    mesh = d.fun.mesh;
    h    = mesh.computeMeanCellSize();

    s.mesh  = mesh;
    s.alpha = 4;
    s.beta  = 0;
    s.theta = 90;
    s.tol0  = 1e-6;
    filter  = NonLinearFilterSegment(s);

    levelSet         = -ones(mesh.nnodes,1);
    sG.backgroundMesh = mesh;
    sG.boundaryMesh   = mesh.createBoundaryMesh();
    uMesh = UnfittedMesh(sG);
    uMesh.compute(levelSet);
    
    ss.mesh    = mesh;
    ss.uMesh   = uMesh;
    ss.filter  = filter;
    ss.epsilon = 3*h;
    ss.value0  = 1;
    pF         = PerimeterFunctional(ss);
    
    P          = pF.computeFunctionAndGradient(d);
end