function testQuadrilateralSolutions
    %PREAMBLE
    % test for the squaredMesh with sideLength = 1 and unitDiv = 3
    sideLength = [1,1];
    theta = [0,90];
    divUnit = 3;
    % Initial Data
    nV = load('nvertQuad.mat');
    bN = load('boundNodesQuad.mat');
    tN = load('totalNodesQuad.mat');
    vC = load('vertCoordQuad.mat');
    b = load('boundaryQuad.mat');
    c = load('coordQuad.mat');
    initialData.c = sideLength;
    initialData.theta = theta;
    initialData.divUnit = divUnit;
    initialData.div = divUnit*sideLength;
    initialData.nvert = nV.nvert;
    initialData.nodes.vert = nV.nvert;
    initialData.nodes.bound = bN.boundNodes;
    initialData.nodes.total = tN.totalNodes;
    initialData.vertCoord = vC.vertCoord;
    initialData.boundCoord = b.boundary;
    initialData.coord = c.coord;
    initialData.filename = 'TestQuadrilateralSol';
    
    %TESTERS
    testers = {'NodesCalculatorTester','VertexCoordinatesCalculatorTester',...
    'BoundaryCoordinatesCalculatorTester','IntersectionCoordComputerTester',...
    'MasterSlaveComputerTester','QuadrilateralNodesCalculatorTester',...
    'NodeCoordinatesComputerTester','MeshCreatorTester'};
    
    for iTest = 1:length(testers)
        Tester.create(testers{iTest},initialData);
    end

end