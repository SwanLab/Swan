function P = computeReferenceGlobalPerimeterEllipse(matFile)
    load(matFile,'d');
    mesh = d.fun.mesh;
    h    = mesh.computeMeanCellSize();

    CAnisotropic = [tand(85), 0; 0, 1/tand(85)];
    aniAlphaDeg = 90;
    R = [cosd(aniAlphaDeg),-sind(aniAlphaDeg)
    sind(aniAlphaDeg), cosd(aniAlphaDeg)];
    CGlobal = R*CAnisotropic*R';

    s.filterType   = 'PDE';
    s.mesh         = mesh;
    s.boundaryType = 'Neumann';
    s.metric       = 'Anisotropy';
    s.trial        = LagrangianFunction.create(mesh,1,'P1');
    s.A            = ConstantFunction.create(CGlobal,mesh);
    filter         = Filter.create(s);
    
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