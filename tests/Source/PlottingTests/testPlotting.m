classdef testPlotting < testNotShowingError...
        & testUnfitted...
        & testLoadStoredVariable
    
    properties (Access = protected, Abstract)
        testName
        meshType
    end
    
    properties (Access = protected)
        variablesToStore = {'coord','connec'};
    end
    
    methods (Access = protected)
        function obj = testPlotting()
            obj.createTopOpt();
            obj.createMesh();
            obj.plot();
        end
        
        function plot(obj)
            if isequal(obj.meshType,'BOUNDARY')            
                obj.mesh.plotBoundary();
            elseif isequal(obj.meshType,'INTERIOR')            
                obj.mesh.plot();
            end
             view(obj.getViewAngle());            
        end
        
        function hasPassed = hasPassed(obj)
            hasPassed = true;
        end
    end
    
    methods (Access = private)
        function angle = getViewAngle(obj)
            if isprop(obj,'viewAngle')
                angle = obj.viewAngle;
            else
                angle = [0 0 1];
            end
        end
    end
end

