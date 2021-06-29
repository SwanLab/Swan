classdef testTriangleToyUntittedExample < PlottingToyUnfittedExample
    
    properties (Access = public)
       testName = 'test_triangleToyUntittedExample';
    end

    methods (Access = public)
        
        function obj = testTriangleToyUntittedExample()
            figure()
            obj.init();
            obj.compute();
            obj.plot();
        end        
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.coord  = [0 0;1 0;1 1;0 1;2 0;2 1;0 2;1 2;2 2; 0.5 0.5; 1.5 0.5; 0.5 1.5; 1.5 1.5];
            obj.connec = [1 2 10; 2 3 10; 10 3 4; 10 4 1; 2 11 3; 2 5 11; 5 6 11; 11 6 3; 3 8 12; 4 3 12; 12 8 7; 12 7 4; 3 6 13; 6 9 13; 13 9 8; 3 13 8];
            obj.levelSet = [-0.05 0.2 -0.5 0.1 0.1 -1 1 -0.2 -0.5 -0.05 -0.05 0.05 -0.5]';
        end             
        
    end
    
end