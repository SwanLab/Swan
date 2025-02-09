classdef NodeCoordinatesComputerTester < Tester
        
    properties (Access = private)
        data
    end
    
    properties (Access = public) 
        testName
        corrValues
        calcValues
    end    
    
    methods (Access = public)
        
        function obj = NodeCoordinatesComputerTester(cParams)
            obj.data = cParams;
            obj.testName = 'NodeCoordinatesComputer';
            obj.loadCorrectValues();
            obj.obtainCalculatedData();
        end
        
    end
    
    methods (Access = protected)
        
        function loadCorrectValues(obj)
            switch obj.data.nvert
                case 4
                    vC = load('vertCoordQuad.mat');
                    bC = load('boundCoordQuad.mat');
                    tC = load('totalCoordQuad.mat');
                case 6
                    vC = load('vertCoordHex.mat');
                    bC = load('boundCoordHex.mat');
                    tC = load('totalCoordHex.mat');
            end
            obj.corrValues(1).Matrix = vC.vertCoord;
            obj.corrValues(2).Matrix = bC.boundary;
            obj.corrValues(3).Matrix = tC.coord;
        end
        
        function obtainCalculatedData(obj)
            solution = NodeCoordinatesComputer(obj.data);
            solution.computeCoordinates();
            obj.calcValues(1).Matrix = solution.vertCoord;
            obj.calcValues(2).Matrix = solution.boundCoord;
            obj.calcValues(3).Matrix = solution.totalCoord;
        end
    
    end
    
end