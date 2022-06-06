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
            switch obj.data.nvert
                case 4
                    vC = load('vertCoordQuad.mat');
                case 6
                    vC = load('vertCoordHex.mat');
            end
            obj.corrValues(1).Matrix = vC.vertCoord;
        end
        
        function obtainCalculatedData(obj)
            solution = VertexCoordinatesCalculator(obj.data);
            obj.calcValues(1).Matrix = solution.vertCoord;
        end
        
    end
    
end