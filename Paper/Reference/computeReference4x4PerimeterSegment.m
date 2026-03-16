function P = computeReference4x4PerimeterSegment(matFile)
    load(matFile,'d');
    mesh = d.fun.mesh;
    h    = mesh.computeMeanCellSize();

    s.mesh  = mesh;
    s.alpha = 4;
    s.beta  = 0;
    s.theta = 90;
    s.tol0  = 1e-6;
    filter  = NonLinearFilterSegment(s);

    xMax = max(mesh.coord(:,1));
    yMax = max(mesh.coord(:,2));
    xMin = min(mesh.coord(:,1));
    yMin = min(mesh.coord(:,2));

    xSide = (xMax-xMin)/4;
    ySide = (yMax-yMin)/4;

    xCG = xMin+xSide/2:xSide:xMax;
    yCG = yMax-ySide/2:-ySide:yMin;

    x0 = repmat(xCG,[1,4]);
    y0 = [repmat(yCG(1),[1,4]),repmat(yCG(2),[1,4]),repmat(yCG(3),[1,4]),repmat(yCG(4),[1,4])];
    
    ss.mesh    = mesh;
    ss.filter  = filter;
    ss.epsilon = 3*h;
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