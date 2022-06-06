classdef VertexCoordinatesCalculatorTester < Tester
    
    properties (Access = private)
        data
    end
    
    properties (Access = protected)
        testName
        corrValues
        calcValues
    end
    
    methods (Access = public)
        
        function obj = VertexCoordinatesCalculatorTester(cParams)
            obj.data = cParams;
            obj.testName = 'VertexCoordinatesCalculatorTester';
            obj.loadCorrectValues();
            obj.obtainCalculatedData();
            obj.verify();
        end
        
    end
    
    methods (Access = private)
        
        function loadCorrectValues(obj)
            vC = load('vertCoordQuad.mat');
            obj.corrValues(1).Matrix = vC.vertCoord;
        end
        
        function obtainCalculatedData(obj)
            solution = VertexCoordinatesCalculator(obj.data);
            obj.calcValues(1).Matrix = solution.vertCoord;
        end
        
    end
    
end