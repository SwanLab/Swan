classdef DiagonalCoordComputerTester < Tester
       
    properties (Access = private)
        data
    end
    
    properties (Access = protected)
        testName
        corrValues
        calcValues
    end
    
    methods (Access = public)
        
        function obj = DiagonalCoordComputerTester(cParams)
            obj.data = cParams;
            obj.testName = 'DiagonalCoordComputer';
            obj.loadCorrectValues();
            obj.obtainCalculatedData();
            obj.verify();
        end
        
    end
    
    methods (Access = protected)
        
        function loadCorrectValues(obj)
            tC = load('totalCoordHex.mat');
            obj.corrValues(1).Matrix = tC.coord;
        end
        
        function obtainCalculatedData(obj)
            solution = DiagonalCoordComputer(obj.data);
            obj.calcValues(1).Matrix = solution.totalCoord;
        end
        
    end
    
end