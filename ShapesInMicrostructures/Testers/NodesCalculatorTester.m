classdef NodesCalculatorTester < Tester
    
    properties (Access = private)
        data
    end
    
    properties (Access = protected) 
        testName
        corrValues
        calcValues
    end    
    
    methods (Access = public)
        
        function obj = NodesCalculatorTester(cParams)
            obj.data = cParams;
            obj.testName = 'NodesCalculatorTester';
            obj.loadCorrectValues();
            obj.obtainCalculatedData();
            obj.verify();
        end
        
    end
    
    methods (Access = protected)
        
        function loadCorrectValues(obj)
            nV = load('nvert.mat');
            bN = load('boundNodes.mat');
            tN = load('totalNodes.mat');
            obj.corrValues(1).Matrix = nV.nsides;
            obj.corrValues(2).Matrix = bN.boundNodes;
            obj.corrValues(3).Matrix = tN.totalNodes;
        end
        
        function obtainCalculatedData(obj)
            solution = NodesCalculator.create(obj.data);
            obj.calcValues(1).Matrix = solution.nvert;
            obj.calcValues(2).Matrix = solution.boundNodes;
            obj.calcValues(3).Matrix = solution.totalNodes;
        end
    
    end
 
end