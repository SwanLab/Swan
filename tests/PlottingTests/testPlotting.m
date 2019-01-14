classdef testPlotting < testNotShowingError...
                        & testUnfitted
%         & testLoadStoredVariable
    
    properties (Access = protected, Abstract)
        testName
    end
    
    methods (Access = protected)
        function obj = testPlotting()
            obj.createTopOpt()
            obj.createMesh()
            obj.plot()
        end
        
        function plot(obj)
            obj.mesh.plot();
        end
        
        function hasPassed = hasPassed(obj)
           hasPassed = false; 
        end
    end
    
end

