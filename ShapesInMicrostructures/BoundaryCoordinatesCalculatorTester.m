classdef BoundaryCoordinatesCalculatorTester < Tester
    
    properties (Access = private)
        data
    end
    
    properties (Access = protected)
        testName
        corrValues
        calcValues
    end
    
    methods (Access = public)
        
        function obj = BoundaryCoordinatesCalculatorTester(cParams)
            obj.data = cParams;
            obj.testName = 'BoundaryCoordinatesCalculatorTester';
            obj.loadCorrectValues();
            obj.obtainCalculatedData();
            obj.verify();
        end
        
    end
    
    methods (Access = private)
        
        function loadCorrectValues(obj)
            bC = load('boundary.mat');
            obj.corrValues(1).Matrix = bC.boundary;
        end
        
        function obtainCalculatedData(obj)
            solution = BoundaryCoordinatesCalculator(obj.data);
            obj.calcValues(1).Matrix = solution.boundCoord;
        end
        
    end
    
end