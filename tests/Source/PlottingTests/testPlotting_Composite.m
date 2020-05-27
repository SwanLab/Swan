classdef testPlotting_Composite < testPlotting
    methods (Access = protected)
        
        function hasPassed = hasPassed(obj)
             coordsOk = isequal(obj.mesh.coord,obj.storedVar{1});
            connecsOk = isequal(obj.mesh.connec,obj.storedVar{2});
            hasPassed = coordsOk && connecsOk;
        end
        
    end
    
end

