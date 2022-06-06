classdef IntersectionCoordComputerTester < Tester
        
    properties (Access = private)
        data
    end
    
    properties (Access = protected) 
        testName
        corrValues
        calcValues
    end    
    
    methods (Access = public)
        
        function obj = IntersectionCoordComputerTester(cParams)
            obj.data = cParams;
            obj.testName = 'IntersectionCoordComputerTester';
            obj.loadCorrectValues();
            obj.obtainCalculatedData();
            obj.verify();
        end
        
    end
    
    methods (Access = protected)
        
        function loadCorrectValues(obj)
            c = load('coordQuad.mat');
            obj.corrValues(1).Matrix = c.coord;
        end
        
        function obtainCalculatedData(obj)
            solution = IntersectionCoordComputer(obj.data);
            obj.calcValues(1).Matrix = solution.totalCoord;
        end
    
    end
    
end