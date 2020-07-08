classdef testPlotting_Composite < testPlotting
    methods (Access = protected)
        
        function hasPassed = hasPassed(obj)
            coordsOk = isequal(obj.mesh,obj.storedVar{1});
            connecsOk = isequal(obj.mesh,obj.storedVar{2});
            hasPassed = coordsOk && connecsOk;
        end
        
    end
    
end

