classdef testPlotting < testNotShowingError...
        & testUnfitted...
        & testLoadStoredVariable
    
    properties (Access = protected, Abstract)
        testName
    end
    
    properties (Access = protected)
        variablesToStore = {'coord','connec'};
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
            if all(size(obj.mesh.coord) == size(obj.storedVar{1})) &&  all(size(obj.mesh.connec) == size(obj.storedVar{2}))
                if all(all(obj.mesh.coord == obj.storedVar{1})) && all(all(obj.mesh.connec == obj.storedVar{2}))
                    hasPassed = true;
                else
                    hasPassed = false;
                    return
                end
            else
                hasPassed = false;
                return
            end
        end
    end
end

