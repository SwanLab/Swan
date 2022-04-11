classdef NodeCoordinatesComputerTester < Tester
        
    properties (Access = private)
        data
    end
    
    properties (Access = protected) 
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
            obj.verify();
        end
        
    end
    
    methods (Access = protected)
        
        function loadCorrectValues(obj)
            vC = load('vertCoord.mat');
            bC = load('boundary.mat');
            tC = load('coord.mat');
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