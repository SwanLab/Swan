classdef testQuadToyUntittedExample < PlottingToyUnfittedExample
    
    properties (Access = public)
        testName = 'test_quadToyUntittedExample';
    end
    
    methods (Access = public)
        
        function obj = testQuadToyUntittedExample()      
            obj.init();
            obj.compute(); 
            obj.plot();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.coord = [0 0;1 0;1 1;0 1;2 0;2 1;0 2;1 2;2 2];
            obj.connec = [1 2 3 4; 2 5 6 3; 4 3 8 7; 3 6 9 8];            
            obj.levelSet = [-0.05 0.2 -0.5 0.1 0.1 -1 1 -0.2 -0.5]';
        end
        
    end
    
end