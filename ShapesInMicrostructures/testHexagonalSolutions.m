function testHexagonalSolutions
    %PREAMBLE
    % test for the hexagonalMesh with sideLength = 1 and unitDiv = 3
    sideLength = [1,1,1];
    theta = [0,60,120];
    divUnit = 3;
    % Initial Data
    nV = load('nvertH.mat');
    bN = load('boundNodesH.mat');
    tN = load('totalNodesH.mat');
    vC = load('vertCoordH.mat');
    b = load('boundCoordH.mat');
    c = load('totalCoordH.mat');
    initialData.c = sideLength;
    initialData.theta = theta;
    initialData.divUnit = divUnit;
    initialData.div = divUnit*sideLength;
    initialData.nvert = nV.nsides;
    initialData.nodes.vert = nV.nsides;
    initialData.nodes.bound = bN.boundNodes;
    initialData.nodes.total = tN.nnodes;
    initialData.vertCoord = vC.vertCoord;
    initialData.boundCoord = b.boundary;
    initialData.totalCoord = c.coord;
    
    testers = {'DiagonalCoordComputerTester','HexagonalNodesCalculatorTester'};
    for iTest = 1:length(testers)
        Tester.create(testers{iTest},initialData);
    end

end