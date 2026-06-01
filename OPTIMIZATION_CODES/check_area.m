clear; clc;
try
    elements = 50;
    fullmesh = TriangleMesh(18,10,elements,elements);
    
    gPar.type          = 'WingShape';
    gPar.xCoorCenter   = -0.01;
    gPar.yCoorCenter   = -0.01;
    gPar.chordRoot     = 7.3;
    gPar.chordTip      = 1.25;
    gPar.semiSpan      = 18.0;
    gPar.sweepDeg      = 25.0;
    g                  = GeometricalFunction(gPar);
    phiFun             = g.computeLevelSetFunction(fullmesh);
    ls                 = phiFun.fValues;

    sUm.backgroundMesh = fullmesh;
    sUm.boundaryMesh   = fullmesh.createBoundaryMesh;
    uMesh              = UnfittedMesh(sUm);
    uMesh.compute(ls);
    wingMesh = uMesh.createInnerMesh();

    area = wingMesh.computeVolume();
    fprintf('COMPUTED_AREA: %.10f\n', area);
catch e
    fprintf('ERROR: %s\n', e.message);
end
