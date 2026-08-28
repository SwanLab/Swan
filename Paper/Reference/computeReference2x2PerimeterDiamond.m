function P = computeReference2x2PerimeterDiamond(degAngle,matFile)
    load(matFile,'d');
    mesh = d.fun.mesh;
    h    = mesh.computeMeanCellSize();
    
    s.mesh        = mesh;
    s.trial       = LagrangianFunction.create(mesh,1,'P1');
    s.senseVector = ConstantFunction.create([0;1],mesh);
    s.ovAngleDeg  = degAngle;
    s.tol         = 3e-4;
    filter        = FilterRegularizedDiamond(s);
    
    xMax = max(mesh.coord(:,1));
    yMax = max(mesh.coord(:,2));
    xMin = min(mesh.coord(:,1));
    yMin = min(mesh.coord(:,2));

    xSide = (xMax-xMin)/2;
    ySide = (yMax-yMin)/2;

    xCG = xMin+xSide/2:xSide:xMax;
    yCG = yMax-ySide/2:-ySide:yMin;

    x0 = repmat(xCG,[1,2]);
    y0 = [repmat(yCG(1),[1,2]),repmat(yCG(2),[1,2])];
    
    ss.mesh    = mesh;
    ss.filter  = filter;
    ss.epsilon = 6*h;
    ss.value0  = 1;
    for i = 1:length(x0)
        ss.uMesh = createBaseDomainPerimeter(mesh,x0(i),y0(i),xSide,ySide);
        pF       = PerimeterFunctional(ss);
        P(i)     = pF.computeFunctionAndGradient(d);
    end
end

function uMesh = createBaseDomainPerimeter(mesh,x0,y0,xSide,ySide)
    s.type             = 'Rectangle';
    s.xCoorCenter      = x0;
    s.yCoorCenter      = y0;
    s.xSide            = xSide;
    s.ySide            = ySide;
    g                  = GeometricalFunction(s);
    lsFun              = g.computeLevelSetFunction(mesh);
    sUm.backgroundMesh = mesh;
    sUm.boundaryMesh   = mesh.createBoundaryMesh();
    uMesh              = UnfittedMesh(sUm);
    uMesh.compute(lsFun.fValues);
end