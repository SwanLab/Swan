function testHexagonalSolutions
    %PREAMBLE
    % test for the hexagonalMesh with sideLength = 1 and unitDiv = 3
    sideLength = [1,1,1];
    theta = [0,60,120];
    divUnit = 3;
    % Initial Data
    nV = load('nvertHex.mat');
    bN = load('boundNodesHex.mat');
    tN = load('totalNodesHex.mat');
    vC = load('vertCoordHex.mat');
    b = load('boundCoordHex.mat');
    c = load('totalCoordHex.mat');
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
    initialData.filename = 'TestHexagonalSol';
    
    %SPECIFIC TESTERS
    testers = {'DiagonalCoordComputerTester','HexagonalNodesCalculatorTester',...
        'MasterSlaveComputerTester'
        };
    
    for iTest = 1:length(testers)
        Tester.create(testers{iTest},initialData);
    end

end