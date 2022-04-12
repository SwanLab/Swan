classdef HexagonalNodesCalculatorTester < Tester
       
    properties (Access = private)
        data
    end
    
    properties (Access = protected)
        testName
        corrValues
        calcValues
    end
    
    methods (Access = public)
        
        function obj = HexagonalNodesCalculatorTester(cParams)
            obj.data = cParams;
            obj.testName = 'HexagonalNodesCalculator';
            obj.loadCorrectValues();
            obj.obtainCalculatedData();
            obj.verify();
        end
        
    end
    
    methods (Access = protected)
        
        function loadCorrectValues(obj)
            bN = load('boundNodesH.mat');
            tN = load('totalNodesH.mat');
            obj.corrValues(1).Matrix = bN.boundNodes;
            obj.corrValues(2).Matrix = tN.nnodes;
        end
        
        function obtainCalculatedData(obj)
            solution = HexagonalNodesCalculator(obj.data);
            obj.calcValues(1).Matrix = solution.boundNodes;
            obj.calcValues(2).Matrix = solution.totalNodes;
        end
        
    end

end